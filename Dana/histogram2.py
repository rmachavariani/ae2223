import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from Dana.class_def import *
from Dana.gridConstruction import getRectangularGridWithVectors
from Dana.gridConstruction import loadData
from Dana.gridConstruction import createVectorObjects
from matplotlib.cm import ScalarMappable
from astropy.convolution import convolve, Gaussian2DKernel
from scipy import stats
import seaborn as sns


def getGeneralDensity():

    data = loadData()
    vectors = createVectorObjects(data)
    vectors = vectors[::1000]

    u_values = []
    v_values = []
    w_values = []

    x_values = []
    y_values = []
    z_values = []

    for vector in vectors:
        u_values.append(vector.u)
        v_values.append(vector.v)
        w_values.append(vector.w)

        x_values.append(vector.x)
        y_values.append(vector.y)
        z_values.append(vector.z)

    # heatmap, xedges, yedges = np.histogram2d(x_values, z_values, bins=40)
    # extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    # plt.imshow(convolve(heatmap, Gaussian2DKernel(stddev=2)), interpolation='none')
    # plt.colorbar(u_values)
    # plt.title('helo')
    # plt.xlabel("X")
    # plt.ylabel("Z")
    # plt.show()

    f = plt.figure()

    plt.rcParams.update({'font.size': 6})

    f, axes = plt.subplots(nrows=2, ncols=1, figsize=(2,3), gridspec_kw={'height_ratios': [1, 2]})
    x_plot = axes[1].scatter(y_values, x_values, c=u_values, edgecolors='none', marker=".",  cmap=cm.jet)

    axes[1].set_xlim([-150,150])
    axes[1].set_ylim([0,250])
    axes[1].yaxis.set_label_position("right")
    axes[1].set_xlabel('y [mm]')
    axes[1].set_ylabel('x [mm]')

    z_plot = axes[0].scatter(y_values, z_values, c=u_values, edgecolors='none', marker=".",  cmap=cm.jet)
    axes[0].set_xlim([-150,150])
    axes[0].set_ylim([0,250])
    axes[0].set_xticks([])
    axes[0].set_ylabel('z [mm]')
    axes[0].yaxis.set_label_position("right")


    bar = f.colorbar(x_plot, ax=axes.ravel().tolist(), location="top")
    bar.ax.set_xlabel('u [m/s]')
    plt.show()


    f, axes = plt.subplots(nrows=2, ncols=1, figsize=(2,3), gridspec_kw={'height_ratios': [1, 2]})
    x_plot = axes[1].scatter(y_values, x_values, c=v_values, edgecolors='none', marker=".",  cmap=cm.jet)
    axes[1].set_xlim([-150,150])
    axes[1].set_ylim([0,250])
    axes[1].set_xlabel('y [mm]')
    axes[1].set_ylabel('x [mm]')
    axes[1].yaxis.set_label_position("right")

    z_plot = axes[0].scatter(y_values, z_values, c=v_values, edgecolors='none', marker=".",  cmap=cm.jet)
    axes[0].set_xlim([-150,150])
    axes[0].set_ylim([0,250])
    axes[0].set_xticks([])
    axes[0].set_ylabel('z [mm]')
    axes[0].yaxis.set_label_position("right")

    bar = f.colorbar(x_plot, ax=axes.ravel().tolist(), location="top")
    bar.ax.set_xlabel('v [m/s]')
    plt.show()

    f, axes = plt.subplots(nrows=2, ncols=1, figsize=(2,3), gridspec_kw={'height_ratios': [1, 2]})
    x_plot = axes[1].scatter(y_values, x_values, c=w_values, edgecolors='none', marker=".",  cmap=cm.jet)
    axes[1].set_xlim([-150,150])
    axes[1].set_ylim([0,250])
    axes[1].set_xlabel('y [mm]')
    axes[1].set_ylabel('x [mm]')
    axes[1].yaxis.set_label_position("right")

    z_plot = axes[0].scatter(y_values, z_values, c=w_values, edgecolors='none', marker=".",  cmap=cm.jet)
    axes[0].set_xlim([-150,150])
    axes[0].set_ylim([0,250])
    axes[0].set_xticks([])
    axes[0].set_ylabel('z [mm]')
    axes[0].yaxis.set_label_position("right")

    bar = f.colorbar(x_plot, ax=axes.ravel().tolist(), location="top")
    bar.ax.set_xlabel('w [m/s]')
    plt.show()





# n_bins = 4
# data = getRectangularGridWithVectors(5000, 5000, 5000)
# data = getSphericalGridWithVectors(50,50)
# bins = getMostPopulated(data, n_bins)
getGeneralDensity()
# histograms2d = get2dHistograms(data,bins)
# histograms = getHistograms(data, bins)