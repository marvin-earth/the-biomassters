from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import imageio as iio
import pandas as pd

from dataset import read_imgs, read_shape_mask

if __name__ == "__main__":
    out_dir = Path("/home/ubuntu/asaph/the-biomassters/1st-place/preds2")
    test_images_dir = Path("/home/ubuntu/asaph/the-biomassters/1st-place/suzano_inputs/suzano_tifs2")
    test_df_path = Path("/home/ubuntu/asaph/the-biomassters/1st-place/suzano_inputs/suzano_meta.csv")

    df = pd.read_csv(test_df_path)
    if "split" in df.columns:
        test_df = df[df.split == "test"].copy()
    else:
        test_df = df
    test_df = test_df.groupby("chip_id").agg(list).reset_index()

    preds = []
    gts = []
    for p in out_dir.glob("*.tif"):
        chip_id = p.name[: len("T8AO03_2023-02-09")]
        fig, ax = plt.subplots(1, 2)
        mask = read_shape_mask(chip_id, test_images_dir)
        res = iio.imread(p)

        ax[0].imshow(res)
        ax[1].imshow(mask)

        res_nan = np.where(mask, res, np.nan)
        pred_agbd = np.nanmean(res_nan)
        pred_agbd_std = np.nanstd(res_nan)
        preds.append(pred_agbd)

        gt = test_df[test_df.chip_id == chip_id].agbd.iloc[0][0]
        gts.append(gt)

        fig.suptitle(f"{chip_id} pred: {pred_agbd:.2f} +- {pred_agbd_std:.2f}    GT: {gt / 0.47:.2f}")
        plt.show()

    preds = np.array(preds)
    gts = np.array(gts)

    preds_filt = preds[preds > 1]
    gts_filt = gts[preds > 1]

    plt.scatter(gts_filt, preds_filt)
    plt.show()
    print("bye")
