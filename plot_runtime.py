import os, sys, math, random, time, csv, copy, argparse

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt

import numpy as np
import cv2

fig = None
ax1 = None
ax2 = None
ax3 = None

def plot_overall():
    # max_vals = [
    #         94.03,
    #         98.06,
    #         94.43,
    #     ]
    # vals = [
    #         [3.6, 4.53, 10.01, 53.98, 14.16],
    #         [4.05, 13.29, 19.0, 27.93, 26.31],
    #         [4.57, 5.05, 11.04, 50.25, 16.29],
    #     ]
    scene_names = ['Dining Room', 'Table', 'Living Room']
    for i, _ in enumerate(vals):
        vals[i].append(max_vals[i]-sum(vals[i]))
    
    vals = np.array(vals).T
    labels = ['Frustum Culling', 'LTC Integration', 'Lights Preproc', 'Occluder Preproc', 'Set Difference', 'Misc']
    colors = ['tomato', 'sandybrown', 'wheat', 'khaki', 'lightgreen', 'lightgrey']

    ind = np.arange(len(max_vals))

    ax1.set_xticks(ind)
    ax1.set_xticklabels(scene_names)
    
    for i in range(vals.shape[0]):
        if i == 0:
            ax1.bar(ind, vals[i, :], 0.35, color=colors[i])
        else:
            ax1.bar(ind, vals[i, :], 0.35, bottom=np.sum(vals[:i, :], axis=0), color=colors[i])

        ax1.set_ylabel('Relative Time')
        ax1.set_xlabel('Scenes')
        ax1.legend(labels=labels)

def plot_lights():
    # max_vals = [
    #         10.01,
    #         19.0,
    #         11.04,
    #     ]
    # vals = [
    #         [0.27, 0.24, 2.66, 0.89, 3.68],
    #         [0.63, 0.57, 2.09, 6.04, 5.51],
    #         [0.42, 0.21, 1.81, 3.28, 2.97],
    #     ]
    scene_names = ['Dining Room', 'Table', 'Living Room']
    for i, _ in enumerate(vals):
        vals[i].append(max_vals[i]-sum(vals[i]))
    
    vals = np.array(vals).T
    labels = ['Clip To Horizon', 'Local Shading Transform', 'Project To Plane', 'Silhouette Comp.', 'Sort Edges', 'Misc']
    colors = ['tomato', 'sandybrown', 'wheat', 'khaki', 'lightgreen', 'lightgrey']

    ind = np.arange(len(max_vals))

    ax2.set_xticks(ind)
    ax2.set_xticklabels(scene_names)

    for i in range(vals.shape[0]):
        if i == 0:
            ax2.bar(ind, vals[i, :], 0.35, color=colors[i])
        else:
            ax2.bar(ind, vals[i, :], 0.35, bottom=np.sum(vals[:i, :], axis=0), color=colors[i])

        ax2.set_ylabel('Relative Time')
        ax2.set_xlabel('Light Preprocess')
        ax2.legend(labels=labels)

def plot_occluder():
    # max_vals = [
    #         53.98,
    #         27.93,
    #         50.25,
    #     ]
    # vals = [
    #         [4.04, 1.38, 10.77, 14.96, 15.65],
    #         [1.11, 0.99, 3.56, 9.64, 7.01],
    #         [1.96, 1.27, 9.03, 15.8, 14.5],
    #     ]
    scene_names = ['Dining Room', 'Table', 'Living Room']
    for i, _ in enumerate(vals):
        vals[i].append(max_vals[i]-sum(vals[i]))
    
    vals = np.array(vals).T
    labels = ['Clip To Horizon', 'Local Shading Transform', 'Project To Plane', 'Silhouette Comp.', 'Sort Edges', 'Misc']
    colors = ['tomato', 'sandybrown', 'wheat', 'khaki', 'lightgreen', 'lightgrey']

    ind = np.arange(len(max_vals))

    ax3.set_xticks(ind)
    ax3.set_xticklabels(scene_names)

    for i in range(vals.shape[0]):
        if i == 0:
            ax3.bar(ind, vals[i, :], 0.35, color=colors[i])
        else:
            ax3.bar(ind, vals[i, :], 0.35, bottom=np.sum(vals[:i, :], axis=0), color=colors[i])

        ax3.set_ylabel('Relative Time')
        ax3.set_xlabel('Occluder Preprocess')
        ax3.legend(labels=labels)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    nrows = 1
    ncols = 3

    fig = plt.figure(figsize=(16, 8))

    ax1 = fig.add_subplot(nrows, ncols, 1)
    ax2 = fig.add_subplot(nrows, ncols, 2)
    ax3 = fig.add_subplot(nrows, ncols, 3)

    plot_overall()
    plot_lights()
    plot_occluder()

    plt.show()