import ee
import geemap
import geopandas as gpd
import mrv.crs
import pandas as pd

from client_utils import read_suzano_inventory
from mrv.satellites.sentinel2 import maskCloudsS2
from mrv.satellites.sentinel1 import get_sentinel1


def get_sentinel2(
    aoi: ee.FeatureCollection,
    start_datetime: str,
    end_datetime: str,
    max_cloud_cover: int = 100,
    spectral_indices=None,
    **mask_kwargs,
):
    col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")

    col = col.filterDate(start_datetime, end_datetime).filterBounds(aoi)

    col = col.filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", max_cloud_cover))

    col = maskCloudsS2(col, **mask_kwargs)

    # col = col.scaleAndOffset()

    if spectral_indices:
        col = col.spectralIndices(spectral_indices)

    return col


def inventory2carbon(vtcc, age, density=0.44, co2eq=False):
    tC = vtcc * density * (1 + -0.033928571 * age + 0.432857143) * (1 + -0.020714286 * age + 0.367142857) * 0.47
    if co2eq:
        return tC * 44 / 12
    else:
        return tC


if __name__ == "__main__":
    ee.Initialize(project="marvin-dmrv-dev")
    data_dir = "/tmp/suzano_data"
    # data_dir = "."

    suzano_inv = read_suzano_inventory().reset_index(["unf"])
    idx = pd.IndexSlice
    suzano_inv = suzano_inv.loc[
        idx[:, "2017-10":],
    ]  # at list 6 months since sentinel-2 started collecting

    suzano_ups = gpd.read_file("/home/ubuntu/asaph/suzano_ups.geojson").set_index("up")

    def get_shape(g, which="both"):
        w, s, e, n = g.bounds
        if which == "both":
            return n - s, e - w
        if which == "height":
            return n - s
        if which == "width":
            return e - w

    suzano_ups["width"] = suzano_ups.to_crs(mrv.crs.PSEUDO_MERCATOR).geometry.apply(lambda g: get_shape(g, "width"))
    suzano_ups["height"] = suzano_ups.to_crs(mrv.crs.PSEUDO_MERCATOR).geometry.apply(lambda g: get_shape(g, "height"))
    suzano_ups = suzano_ups.loc[(suzano_ups.width < 2560) & (suzano_ups.height < 2560)]
    suzano_ups["area"] = suzano_ups.to_crs(mrv.crs.PSEUDO_MERCATOR).area / 1e4

    shared_idx = list(set(suzano_inv.index.unique(level="up")) & set(suzano_ups.index.unique(level="up")))
    suzano_ups = suzano_ups.loc[shared_idx].sort_index()
    suzano_inv = suzano_inv.loc[shared_idx].sort_index()

    suzano_inv = suzano_inv[["unf", "measure_date", "age", "vtcc"]]
    suzano_inv.loc[:, "agbd"] = suzano_inv.apply(lambda x: inventory2carbon(x.vtcc, x.age, co2eq=False), axis=1)

    for up in suzano_ups.sort_values("area", ascending=False).head(10).index:
        print(up, suzano_ups["area"].loc[up])
        ee_poly = geemap.gdf_to_ee(suzano_ups.loc[[up]])
        for mdate in suzano_inv.loc[up, "measure_date"]:
            start_date = mdate + pd.Timedelta(-180, "D")
            for i in range(12):
                s2 = get_sentinel2(ee_poly, f"{start_date:%Y-%m-%d}", f'{start_date + pd.Timedelta(30, "D"):%Y-%m-%d}')
                s2_median = s2.select(
                    [
                        "B2",
                        "B3",
                        "B4",
                        "B5",
                        "B6",
                        "B7",
                        "B8",
                        "B8A",
                        "B11",
                        "B12",
                        "MSK_CLDPRB",
                    ]
                ).median()

                # s2_median_np = geemap.ee_to_numpy(
                #     s2_median,
                #     region=ee_poly,
                #     scale=10
                #     # crs=s2.first().select("B2").projection().crs().getInfo(),
                # )
                # print(s2_median_np.shape)

                geemap.download_ee_image(
                    s2_median,
                    f"{data_dir}/{up}_{mdate:%Y-%m-%d}_S2_{i:0>2}.tif",
                    region=ee_poly.geometry(),
                    scale=10,
                    crs=s2.first().select("B2").projection().crs().getInfo(),
                )

                # todo: there are almost no ascending S1 over Brazil (true?), always ignore or add .size().getInfo()?
                # s1_ascending = get_sentinel1(
                #     ee_poly,
                #     f"{start_date:%Y-%m-%d}",
                #     f'{start_date + pd.Timedelta(30, "D"):%Y-%m-%d}',
                #     orbit='ASCENDING',
                # ).median()

                s1_descending = get_sentinel1(
                    ee_poly,
                    f"{start_date:%Y-%m-%d}",
                    f'{start_date + pd.Timedelta(30, "D"):%Y-%m-%d}',
                    orbit="DESCENDING",
                )

                s1_descending_median = s1_descending.median().select(["VV", "VH"], ["VV_DESCENDING", "VH_DESCENDING"])

                # s1 = s1_ascending.addBands(s1_descending)
                if s1_descending.size().getInfo() > 0:
                    geemap.download_ee_image(
                        s1_descending_median,
                        f"{data_dir}/{up}_{mdate:%Y-%m-%d}_S1_{i:0>2}.tif",
                        region=ee_poly.geometry(),
                        scale=10,
                        crs=s1_descending.first().select("VV").projection().crs().getInfo(),
                    )

                start_date += pd.Timedelta(30, "D")
