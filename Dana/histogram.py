import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from Dana.class_def import *
from Dana.gridConstruction import getRectangularGridWithVectors
from Dana.gridConstruction import loadData
from Dana.gridConstruction import createVectorObjects



def getMostPopulated(data, n_bins):

    sizeX = np.size(data,axis=0)
    sizeY = np.size(data,axis=1)
    sizeZ = np.size(data,axis=2)

    lenArray = []

    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                lenArray.append([i, j, k, len(data[i][j][k].vectors)])


    lenArray = np.array(lenArray)
    largestKeys = np.argpartition(lenArray[:,3], -n_bins)[-n_bins:]
    bins = lenArray[largestKeys]

    return bins

def getGeneralDensity():

    data = loadData()
    vectors = createVectorObjects(data)

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

    plt.hexbin(y_values, x_values, C=u_values, gridsize=50, cmap=cm.jet, bins=None)
    plt.axis([ min(y_values), max(y_values), min(x_values), max(x_values),])
    plt.xlabel('y[mm]')
    plt.ylabel('x[mm]')
    cb = plt.colorbar()
    cb.set_label('u [m/s]')
    plt.title("Velocity Heatmap in x-y for u")
    plt.show()

    plt.hexbin(y_values, z_values, C=u_values, gridsize=50, cmap=cm.jet, bins=None)
    plt.axis([ min(y_values), max(y_values), min(z_values),max(z_values)])
    plt.xlabel('y[mm]')
    plt.ylabel('z[mm]')
    cb = plt.colorbar()
    cb.set_label('u [m/s]')
    plt.title("Velocity Heatmap in z-y for u")
    plt.show()

    plt.hexbin(y_values, x_values, C=v_values, gridsize=50, cmap=cm.jet, bins=None)
    plt.axis([ min(y_values), max(y_values), min(x_values), max(x_values),])
    plt.xlabel('y[mm]')
    plt.ylabel('x[mm]')
    cb = plt.colorbar()
    cb.set_label('v [m/s]')
    plt.title("Velocity Heatmap in x-y for v")
    plt.show()

    plt.hexbin(y_values, z_values, C=v_values, gridsize=50, cmap=cm.jet, bins=None)
    plt.axis([ min(y_values), max(y_values), min(z_values),max(z_values)])
    plt.xlabel('y[mm]')
    plt.ylabel('z[mm]')
    cb = plt.colorbar()
    cb.set_label('v [m/s]')
    plt.title("Velocity Heatmap in z-y for v")
    plt.show()

    plt.hexbin(y_values, x_values, C=w_values, gridsize=50, cmap=cm.jet, bins=None)
    plt.axis([ min(y_values), max(y_values), min(x_values), max(x_values),])
    plt.xlabel('y[mm]')
    plt.ylabel('x[mm]')
    cb = plt.colorbar()
    cb.set_label('w [m/s]')
    plt.title("Velocity Heatmap in x-y for w")
    plt.show()

    plt.hexbin(y_values, z_values, C=w_values, gridsize=50, cmap=cm.jet, bins=None)
    plt.axis([ min(y_values), max(y_values), min(z_values),max(z_values)])
    plt.xlabel('y[mm]')
    plt.ylabel('z[mm]')
    cb = plt.colorbar()
    cb.set_label('w [m/s]')
    plt.title("Velocity Heatmap in z-y for w")
    plt.show()





def getHistograms(data, bins):

    for bin in bins:
        bin = data[bin[0], bin[1], bin[2]]
        vectors = bin.vectors
        u_values = []
        v_values = []
        w_values = []

        for vector in vectors:
            u_values.append(vector.u)
            v_values.append(vector.v)
            w_values.append(vector.w)


        plt.hist(u_values, bins='auto', alpha=0.5, label="U",  rwidth=0.85, histtype='step', fc='none', lw=1.5)
        plt.hist(v_values,  bins='auto', alpha=0.5, label="V",  rwidth=0.85, histtype='step', fc='none', lw=1.5)
        plt.hist(w_values,  bins='auto', alpha=0.5, label="W", rwidth=0.85, histtype='step', fc='none', lw=1.5)
        plt.legend(loc='upper left')


        plt.title("Histogram with 'auto' bins")
        plt.show()


def get2dHistograms(data, bins):
    for bin in bins:
        bin = data[bin[0], bin[1], bin[2]]
        vectors = bin.vectors
        u_values = []
        v_values = []
        w_values = []

        x_values = []
        y_values = []

        for vector in vectors:
            u_values.append(vector.u)
            v_values.append(vector.v)
            w_values.append(vector.w)

            x_values.append(vector.x)
            y_values.append(vector.y)

        plt.hexbin(x_values, y_values, C=u_values, gridsize=30, cmap=cm.jet, bins=None)
        plt.axis([min(x_values), max(x_values), min(y_values), max(y_values)])

        cb = plt.colorbar()
        cb.set_label('mean value')

        plt.show()


n_bins = 10
data = getRectangularGridWithVectors(5000, 5000, 5000)
bins = getMostPopulated(data, n_bins)
# getGeneralDensity()
# histograms2d = get2dHistograms(data,bins)
histograms = getHistograms(data, bins)