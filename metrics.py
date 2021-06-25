import os, sys, math, random, time, csv, copy, argparse

from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma)

import numpy as np
import cv2

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--img_1', type=str, default='')
    parser.add_argument('--img_2', type=str, default='')

    args = parser.parse_args()

    im1 = cv2.imread('%s' % args.img_1).astype(np.float) / 255.0
    im2 = cv2.imread('%s' % args.img_2).astype(np.float) / 255.0

    error = np.abs(im1 - im2)
    mae = np.mean(error)

    error = (im1 - im2) ** 2
    mse = np.mean(error)
    rmse = np.sqrt(mse)

    print("MAE: %f" % mae)
    print("RMSE: %f" % rmse)