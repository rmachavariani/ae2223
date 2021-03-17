import numpy as np
from class_def import *
from math import floor, sqrt
import time


def loadData():
    ''' Loads data from carMirrorData.dat
    :return: nrVectors x 6 np array
    '''
    # load the data
    t1 = time.time()
    data = np.loadtxt("carMirrorData.dat", max_rows = 500)
    t2 = time.time()
    print("Loading done in ", "{:.2f}".format(t2 - t1), " s")

    return data


def determineMaxMin(data):
    '''Determines the mininum and maximum value for every dimension
    :return: xMin, xMax, yMin, yMax, zMin, zMax
    '''

    t1 = time.time()

    # determine min and max
    xMin = np.amin(data[:, 0])
    xMax = np.amax(data[:, 0])
    yMin = np.amin(data[:, 1])
    yMax = np.amax(data[:, 1])
    zMin = np.amin(data[:, 2])
    zMax = np.amax(data[:, 2])

    # report to user
    t2 = time.time()
    print("Max and min found in ", "{:.2f}".format(t2 - t1), " s")

    return xMin, xMax, yMin, yMax, zMin, zMax


def createVectorObjects(data):
    ''' Creates objects from the particle/vector data
    :param data: raw data in numpy array
    :return: 1D numpy array with vectors as vector objects
    '''

    t1 = time.time()

    # create empty numpy array
    dataPoints = np.empty(np.size(data, axis=0), dtype=object)

    # loop over data and create vector object for each particle row
    for i in range(np.size(data, axis=0)):
        dataPoints[i] = vector(data[i, :])

    # report to user
    t2 = time.time()
    print("Objects created in ", "{:.2f}".format(t2 - t1), " s")

    return dataPoints


def createGridPitchAndRadius(pitch, radius, xMin, xMax, yMin, yMax, zMin, zMax):

    t1 = time.time()

    # calculate amount of bins in every direction
    nrBinsX = floor((xMax - xMin) / pitch) + 2
    nrBinsY = floor((yMax - yMin) / pitch) + 2
    nrBinsZ = floor((zMax - zMin) / pitch) + 2


    # create empty 3D array
    grid = np.empty((nrBinsX, nrBinsY, nrBinsZ), dtype=object)

    # set radius of bins and
    gridBin.radius = radius
    gridBin.nrBinsX = nrBinsX
    gridBin.nrBinsY = nrBinsY
    gridBin.nrBinsZ = nrBinsZ

    # define x, y and z coordinates of center bin
    x = np.array([(xMin + i * pitch) for i in range(nrBinsX)])
    y = np.array([(yMin + i * pitch) for i in range(nrBinsY)])
    z = np.array([(zMin + i * pitch) for i in range(nrBinsZ)])

    # fill matrix with bin objects by looping over matrix
    for i in range(nrBinsX):
        for j in range(nrBinsY):
            for k in range(nrBinsZ):
                grid[i, j, k] = gridBin(x[i], y[j], z[k])

    # report to user
    t2 = time.time()
    print('Grid created in ', "{:.2f}".format(t2 - t1), " s")

    return grid

def createGrid(nrBinsX, nrBinsY, nrBinsZ, radius, xMin, xMax, yMin, yMax, zMin, zMax):
    '''Creates the grid by generating a 3D numpy array filled with
    objects of class gridbin
    :return: nrX x nrY x nrZ numpy array
    '''

    t1 = time.time()

    # create empty 3D array
    grid = np.empty((nrBinsX, nrBinsY, nrBinsZ), dtype=object)

    # calculate width of bins in all directions
    widthX = (xMax - xMin) / nrBinsX
    widthY = (yMax - yMin) / nrBinsY
    widthZ = (zMax - zMin) / nrBinsZ

    # set widths of bin to class static members
    gridBin.widthX = widthX
    gridBin.widthY = widthX
    gridBin.widthZ = widthZ

    # define x, y and z coordinates of center bin
    x = np.linspace(xMin, xMax - widthX, nrBinsX) + widthX / 2
    y = np.linspace(yMin, yMax - widthY, nrBinsY) + widthY / 2
    z = np.linspace(zMin, zMax - widthZ, nrBinsZ) + widthZ / 2

    # fill matrix with bin objects by looping over matrix
    for i in range(nrBinsX):
        for j in range(nrBinsY):
            for k in range(nrBinsZ):
                grid[i, j, k] = gridBin(x[i], y[j], z[k])

    # report to user
    t2 = time.time()
    print('Grid created in ', "{:.2f}".format(t2 - t1), " s")

    return grid


def assignVectorsToGrid(vectors, grid):

    t1 = time.time()

    # get bin radius and amount of bins in each direction
    radius = gridBin.radius
    nrBinsX = gridBin.nrBinsX
    nrBinsY = gridBin.nrBinsY
    nrBinsZ = gridBin.nrBinsZ

    # loop through all bins
    for i in range(nrBinsX):
        for j in range(nrBinsY):
            for k in range(nrBinsZ):

                # get bin object
                aBin = grid[i][j][k]

                # get coordinates
                x = aBin.x
                y = aBin.y
                z = aBin.z

                # loop over all vectors
                for vector in vectors:

                    # get coordinates
                    xx = vector.x
                    yy = vector.y
                    zz = vector.z

                    if sqrt((xx-x)**2 + (yy-y)**2 + (zz-z)**2) <= radius:

                        aBin.addVector(vector)

    # report to user
    t2 = time.time()
    print("Assigning of vectors to bins completed in ", "{:.2f}".format(t2 - t1), " s")

    return grid
#-------------------------------MAIN--------------------------------#

def getSphericalGridWithVectors(pitch,radius):

    t1 = time.time()

    # load the data
    data = loadData()

    # determine max and min of data in every dimension
    minMax = determineMaxMin(data)

    # set parameters for bins
    xMin = minMax[0]
    xMax = minMax[1]
    yMin = minMax[2]
    yMax = minMax[3]
    zMin = minMax[4]
    zMax = minMax[5]

    # transform raw data into vector objects
    dataPoints = createVectorObjects(data)

    # create bins in grid
    grid = createGridPitchAndRadius(pitch,radius,xMin,xMax,yMin,yMax,zMin,zMax)

    # assign vector objects to correct bins
    # grid is the 3D array filled with gridBin objects containing
    # the correct vector objects
    grid = assignVectorsToGrid(dataPoints,grid)

    # report to user
    t2 = time.time()
    print("Total time: ","{:.2f}".format(t2-t1)," s")

    return grid

grid = getSphericalGridWithVectors(10,10)