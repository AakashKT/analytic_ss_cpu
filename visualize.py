import os, sys, math, random, time, csv, copy, argparse

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt

import numpy as np
import cv2

fig = None
ax1 = None
ax2 = None
data = {}
showOccluders = False
showSpherical = True

prev_mousex = 0
prev_mousey = 0

def plotXY(pixel):
    global ax1, data, showOccluders

    ax1.clear()
    ax2.clear()

    pix_str = str(pixel[0]) + '_' + str(pixel[1])
    pix = copy.deepcopy(data[pix_str])

    num_lights = int(pix[0])
    pix = pix[1:]

    for l in range(num_lights):
        ne = int(pix[0])
        pix = pix[1:]

        phi_theta = pix[:ne*2]
        pix = pix[ne*2:]

        xy = pix[:ne*3]
        pix = pix[ne*3:]

        xy = xy + xy[0:3]
        xy = np.array(xy).reshape(ne+1, 3)

        phi_theta = phi_theta + phi_theta[0:2]
        phi_theta = np.array(phi_theta).reshape(ne+1, 2)

        if showSpherical:
            ax1.set_xlim(0, 2*np.pi, auto=False)
            ax1.set_ylim(0, np.pi, auto=False)
            ax1.invert_yaxis()
            ax1.set_aspect('equal', adjustable='box')

            ax1.plot(phi_theta[:, 0], phi_theta[:, 1], color='b')
            # ax1.plot([0, 2*np.pi], [np.pi/2.0, np.pi/2.0], color='k')
        else:
            ax1.plot(xy[:, 0], xy[:, 1], color='b')

        # Occluders
        num_occ = int(pix[0])
        pix = pix[1:]

        for o in range(num_occ):            
            ne = int(pix[0])
            pix = pix[1:]

            phi_theta = pix[:ne*2]
            pix = pix[ne*2:]

            xy = pix[:ne*3]
            pix = pix[ne*3:]

            xy = xy + xy[0:3]
            xy = np.array(xy).reshape(ne+1, 3)

            phi_theta = phi_theta + phi_theta[0:2]
            phi_theta = np.array(phi_theta).reshape(ne+1, 2)

            if showOccluders:
                if showSpherical:
                    ax1.plot(phi_theta[:, 0], phi_theta[:, 1], color='r')
                    ax2.plot(phi_theta[:, 0], phi_theta[:, 1], color='r', alpha=0.3)
                else:
                    ax1.plot(xy[:, 0], xy[:, 1], color='r')
                    ax2.plot(xy[:, 0], xy[:, 1], color='r', alpha=0.3)
        
        # Clipped
        # print(pix)
        num_clip = int(pix[0])
        pix = pix[1:]

        for o in range(num_clip):
            # print(pix[0])
            ne = int(pix[0])
            pix = pix[1:]

            phi_theta = pix[:ne*2]
            pix = pix[ne*2:]

            xy = pix[:ne*3]
            pix = pix[ne*3:]

            xy = xy + xy[0:3]
            xy = np.array(xy).reshape(ne+1, 3)

            phi_theta = phi_theta + phi_theta[0:2]
            phi_theta = np.array(phi_theta).reshape(ne+1, 2)

            if showSpherical:
                ax2.set_xlim(0, 2*np.pi, auto=False)
                ax2.set_ylim(0, np.pi, auto=False)
                ax2.invert_yaxis()
                ax2.set_aspect('equal', adjustable='box')

                ax2.plot(phi_theta[:, 0], phi_theta[:, 1], color='g')
                # ax1.plot([0, 2*np.pi], [np.pi/2.0, np.pi/2.0], color='k')
            else:
                ax2.plot(xy[:, 0], xy[:, 1], color='g')
    
    fig.canvas.draw()
    ax1.figure.canvas.draw()
    ax2.figure.canvas.draw()

def onclick(event):
    global prev_mousex, prev_mousey

    # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
    #       ('double' if event.dblclick else 'single', event.button,
    #        event.x, event.y, event.xdata, event.ydata))
    
    prev_mousex = int(event.xdata)
    prev_mousey = int(event.ydata)
    
    plotXY(( int(event.xdata), int(event.ydata) ))

def on_key(event):
    global showOccluders, showSpherical

    # print('you pressed', event.key, event.xdata, event.ydata)

    if event.key == 'b' or event.key == 'B':
        showOccluders = not showOccluders
        plotXY((prev_mousex, prev_mousey))
    
    elif event.key == 'h' or event.key == 'H':
        showSpherical = not showSpherical
        plotXY((prev_mousex, prev_mousey))
    
    elif event.key == 'z' or event.key == 'Z':
        extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig('sph_poly.png', bbox_inches=extent, transparent=True, dpi=500)

        extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig('sph_pol_1.png', bbox_inches=extent, transparent=True, dpi=500)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--pixel_file', type=str, default='Pixel data from PBRT')
    parser.add_argument('--image_file', type=str, default='Pixel data from PBRT')

    args = parser.parse_args()

    f = csv.reader(open(args.pixel_file))
    for l in f:
        temp = [float(i) for i in l]
        data[str( int(temp[0]) ) + '_' + str( int(temp[1])) ] = temp[2:]

    nrows = 1
    ncols = 3

    fig = plt.figure(figsize=(20, 15))
    button_press_id = fig.canvas.mpl_connect('button_press_event', onclick)
    key_press_id = fig.canvas.mpl_connect('key_press_event', on_key)

    ax1 = fig.add_subplot(nrows, ncols, 2)
    ax1.set_xlim(0, 2*np.pi, auto=False)
    ax1.set_ylim(0, np.pi, auto=False)
    ax1.invert_yaxis()
    ax1.set_aspect('equal', adjustable='box')

    ax2 = fig.add_subplot(nrows, ncols, 3)
    ax2.set_aspect('equal', adjustable='box')

    # Plot everything
    # plotXY((0, 0))

    render = cv2.imread(args.image_file)

    img_ax = fig.add_subplot(nrows, ncols, 1)
    img_ax.imshow(cv2.cvtColor(render, cv2.COLOR_BGR2RGB))

    plt.show()