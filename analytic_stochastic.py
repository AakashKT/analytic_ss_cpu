import os, sys, math, random, time, csv, copy, argparse

from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma)

import numpy as np
import cv2

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--scene_dir', type=str, default='')

    args = parser.parse_args()

    print('===================================================================')
    print('DO NOT INCLUDE BELOW IN RUN-TIME')
    print('===================================================================')
    os.system('build/pbrt %s/scene_un.pbrt --outfile %s/un.png' % (args.scene_dir, args.scene_dir))
    print('===================================================================')
    print('DO NOT INCLUDE ABOVE IN RUN-TIME')
    print('===================================================================')
    print('')


    os.system('build/pbrt %s/scene_sn.pbrt --outfile %s/sn.png' % (args.scene_dir, args.scene_dir))
    os.system('build/pbrt %s/scene_ltc.pbrt --outfile %s/ltc.png' % (args.scene_dir, args.scene_dir))

    ltc = cv2.imread('%s/ltc.png' % args.scene_dir).astype(np.float) / 255.0
    sn = cv2.imread('%s/sn.png' % args.scene_dir).astype(np.float) / 255.0
    un = cv2.imread('%s/un.png' % args.scene_dir).astype(np.float) / 255.0
    # sn = cv2.imread('%s/sn.png' % args.scene_dir)
    # un = cv2.imread('%s/un.png' % args.scene_dir)

    time_ = time.time()

    sn = denoise_bilateral(sn, multichannel=True)
    un = denoise_bilateral(un, multichannel=True)
    # sn = cv2.bilateralFilter(sn, 4, 75, 75).astype(np.float) / 255.0
    # un = cv2.bilateralFilter(un, 4, 75, 75).astype(np.float) / 255.0

    ratio = np.clip(sn/un, 0, 1)
    output = ltc * ratio

    time_ = time.time() - time_
    print('===================================================================')
    print('TIME FOR DENOISING: %f' % time_)
    print('===================================================================')

    output = (output * 255.0).astype(np.uint8)
    cv2.imwrite('%s/analytic_stochastic.png' % args.scene_dir, output)