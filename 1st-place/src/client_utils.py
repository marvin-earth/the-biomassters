import io

import pandas as pd

from mrv.utils import read_file_from_s3


def read_suzano_inventory(unf="ms"):
    df_inventory = pd.read_excel(
        io.BytesIO(
            read_file_from_s3(
                "marvin-client-data", f"suzano/inventory/{unf.lower()}/bd_inventario_{unf.lower()}_marvin_filtrado.xlsx"
            )
        ),
        parse_dates=["data_plantio", "data_medicao", "data_fim_corte"],
    )
    df_inventory = df_inventory.rename(
        columns={
            "ciclo": "cycle",
            "rotacao": "rotation",
            "n_parcela": "n_parcel",
            "idade_inv": "age",
            "vtcc_viva": "vtcc",
            "fuste_viva": "num_live_trees",
            "imaccc": "productivity",
            "data_plantio": "plant_date",
            "data_medicao": "measure_date",
            "data_fim_corte": "harvest_date",
        }
    )
    df_inventory = df_inventory[
        [
            "unf",
            "up",
            "measure_date",
            "age",
            "vtcc",
            "productivity",
            "num_live_trees",
            "plant_date",
            "harvest_date",
            "matgen",
            "cycle",
            "rotation",
            "lat",
            "lon",
            "n_parcel",
            "macroambiente",
        ]
    ]
    df_inventory["plant_date"] = df_inventory["plant_date"].dt.tz_localize(None).dt.normalize()
    df_inventory["measure_date"] = df_inventory["measure_date"].dt.tz_localize(None).dt.normalize()
    df_inventory["harvest_date"] = df_inventory["harvest_date"].dt.tz_localize(None).dt.normalize()
    df_inventory_means = df_inventory.groupby(["unf", "up", "plant_date", "measure_date"]).mean(numeric_only=True)
    df_inventory_means.loc[:, ["harvest_date", "matgen", "macroambiente"]] = df_inventory.groupby(
        ["unf", "up", "plant_date", "measure_date"]
    )[["harvest_date", "matgen", "macroambiente"]].first()
    df_inventory_means = df_inventory_means.reset_index()

    # remove invalid measure dates
    df_inventory_means["valid_dates"] = (df_inventory_means.measure_date >= df_inventory_means.plant_date) & (
        (df_inventory_means.measure_date <= df_inventory_means.harvest_date) | df_inventory_means.harvest_date.isna()
    )
    # valid_indices = df_inventory_means.groupby(["unf", "up", "plant_date"])["valid_dates"].all()
    df_inventory = df_inventory_means.set_index(["unf", "up", "plant_date"]).loc[
        df_inventory_means.groupby(["unf", "up", "plant_date"])["valid_dates"].all()
    ]

    return df_inventory


def up2unf(up):
    up = up.lower()
    if up[0] in ["c", "t"]:
        return "MS"
    elif up[0] in ["s", "n", "j"]:
        return "SP"
    elif up[0] in ["e", "b", "m"]:
        return "BAES"
    elif up[0] in ["a", "i", "o", "u", "p", "l"]:
        return "MA"
    else:
        raise ValueError(f"unknown parcel code: {up}")


def up2farm(up):
    return up[:4]
