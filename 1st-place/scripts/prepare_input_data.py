from pathlib import Path

import ee
import geemap
import geopandas as gpd
import mrv.crs
import pandas as pd
from shapely import box

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
    data_dir = Path("/home/ubuntu/asaph/the-biomassters/1st-place/suzano_inputs/suzano_tifs2")

    data_dir.mkdir(exist_ok=True, parents=True)
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

    subset = suzano_ups.sort_values("area", ascending=False).head(10)
    subset_inv = suzano_inv.loc[subset.index]
    subset_inv["chip_id"] = (
        subset_inv.index.get_level_values("up") + "_" + subset_inv["measure_date"].dt.strftime("%Y-%m-%d")
    )
    subset_inv.to_csv("/home/ubuntu/asaph/the-biomassters/1st-place/suzano_inputs/suzano_meta.csv")

    bounds = suzano_ups.to_crs(mrv.crs.PSEUDO_MERCATOR).bounds
    bounds.loc[:, ["width"]] = bounds["maxx"] - bounds["minx"]
    bounds.loc[:, ["height"]] = bounds["maxy"] - bounds["miny"]
    bounds.loc[:, ["maxx_padded"]] = bounds["maxx"] + ((2560 - bounds["width"]) / 2)
    bounds.loc[:, ["minx_padded"]] = bounds["minx"] - ((2560 - bounds["width"]) / 2)
    bounds.loc[:, ["maxy_padded"]] = bounds["maxy"] + ((2560 - bounds["height"]) / 2)
    bounds.loc[:, ["miny_padded"]] = bounds["miny"] - ((2560 - bounds["height"]) / 2)
    bounds["boxes"] = bounds.apply(
        lambda r: box(r["minx_padded"], r["miny_padded"], r["maxx_padded"], r["maxy_padded"]), axis=1
    )
    bounds = gpd.GeoDataFrame(bounds, geometry="boxes", crs=mrv.crs.PSEUDO_MERCATOR)
    bounds = bounds.to_crs(mrv.crs.WGS84)

    for up in subset.index:
        print(up, suzano_ups["area"].loc[up])
        ee_poly = geemap.gdf_to_ee(suzano_ups.loc[[up]])
        ee_bbox = ee.Geometry.BBox(*bounds.loc[up].boxes.bounds)
        for mdate in suzano_inv.loc[up, "measure_date"]:
            start_date = mdate + pd.Timedelta(-180, "D")

            for i in range(12):

                # Sentinel-2
                s2 = get_sentinel2(
                    ee_poly,
                    f"{start_date:%Y-%m-%d}",
                    f'{start_date + pd.Timedelta(30, "D"):%Y-%m-%d}',
                    max_cloud_cover=50,
                )

                if i == 0:
                    geemap.download_ee_image(
                        ee_poly.reduceToImage(properties=["area"], reducer=ee.Reducer.first()).gt(0),
                        f"{data_dir}/{up}_{mdate:%Y-%m-%d}_MASK.tif",
                        region=ee_bbox,
                        scale=10,
                        crs=s2.first().select("B2").projection().crs().getInfo(),
                    )

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

                geemap.download_ee_image(
                    s2_median,
                    f"{data_dir}/{up}_{mdate:%Y-%m-%d}_S2_{i:0>2}.tif",
                    region=ee_bbox,
                    scale=10,
                    crs=s2.first().select("B2").projection().crs().getInfo(),
                    # unmask_value=0,  todo: add unmask value?
                )

                # Sentinel-1
                s1_desc = get_sentinel1(
                    ee_poly,
                    f"{start_date:%Y-%m-%d}",
                    f'{start_date + pd.Timedelta(30, "D"):%Y-%m-%d}',
                    orbit="DESCENDING",
                )
                if s1_desc.size().getInfo() > 0:
                    s1_desc_median = s1_desc.median().select(["VV", "VH"], ["VV_DESCENDING", "VH_DESCENDING"])

                    s1_asc = get_sentinel1(
                        ee_poly,
                        f"{start_date:%Y-%m-%d}",
                        f'{start_date + pd.Timedelta(30, "D"):%Y-%m-%d}',
                        orbit="ASCENDING",
                    )
                    if s1_asc.size().getInfo() > 0:
                        s1_asc_median = s1_asc.median().select(["VV", "VH"], ["VV_ASCENDING", "VH_ASCENDING"])
                    else:
                        s1_asc_median = s1_desc_median.select(
                            ["VV_DESCENDING", "VH_DESCENDING"], ["VV_ASCENDING", "VH_ASCENDING"]
                        )

                    s1 = s1_desc_median.addBands(s1_asc_median)

                    geemap.download_ee_image(
                        s1,
                        f"{data_dir}/{up}_{mdate:%Y-%m-%d}_S1_{i:0>2}.tif",
                        region=ee_bbox,
                        scale=10,
                        crs=s1_desc.first().select("VV").projection().crs().getInfo(),
                        unmask_value=-9999,
                    )

                    start_date += pd.Timedelta(30, "D")
