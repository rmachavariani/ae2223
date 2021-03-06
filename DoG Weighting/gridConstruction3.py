import numpy as np
from class_def import *
from math import floor
import time


def loadData():
    ''' Loads data from carMirrorData.dat
    :return: nrVectors x 6 np array
    '''
    # load the data
    t1 = time.time()
    data = np.loadtxt("/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/My coding/carMirrorData.dat", max_rows=20000)
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


def createGrid(nrBinsX, nrBinsY, nrBinsZ, xMin, xMax, yMin, yMax, zMin, zMax):
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

    # set number of bins in class
    gridBin.nrBinsX = nrBinsX
    gridBin.nrBinsY = nrBinsY
    gridBin.nrBinsZ = nrBinsZ

    # set min and max values
    gridBin.xMin = xMin
    gridBin.xMax = xMax
    gridBin.yMin = yMin
    gridBin.yMax = yMax
    gridBin.zMin = zMin
    gridBin.zMax = zMax

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
    ''' Assigns all the vectors to correct entry of the grid
    :param vectors: 1D numpy array with vector objects
    :param grid: 3D numpy array with gridBin objects
    :return: 3D numpy array with gridBin objects with corresponding
    vector objects attached in a list
    '''

    t1 = time.time()

    # get bin widths
    widthX = gridBin.widthX
    widthY = gridBin.widthY
    widthZ = gridBin.widthZ

    # loop through all vectors
    for vector in vectors:

        # get coordinates vector
        x = vector.x
        y = vector.y
        z = vector.z

        # calculate index
        iX = int(floor((x - gridBin.xMin) / widthX))
        iY = int(floor((y - gridBin.yMin) / widthY))
        iZ = int(floor((z - gridBin.zMin) / widthZ))

        # correct for the edge values
        if iX == gridBin.nrBinsX:
            iX += -1

        if iY == gridBin.nrBinsY:
            iY += -1

        if iZ == gridBin.nrBinsZ:
            iZ += -1

        # assign to correct bin
        grid[iX, iY, iZ].addVector(vector)

    # report to user
    t2 = time.time()
    print("Assigning of vectors to bins completed in ", "{:.2f}".format(t2 - t1), " s")

    return grid


def showAmountOfVectorsInBin(grid):
    ''' reports on amount of vectors in every bin'''

    sizeX = np.size(grid, axis=0)
    sizeY = np.size(grid, axis=1)
    sizeZ = np.size(grid, axis=2)

    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                print(len(grid[i][j][k].vectors))


# ------------------------------MAIN-----------------------------#

def getRectangularGridWithVectors(nrBinsX, nrBinsY, nrBinsZ):
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
    print(zMax)

    # transform raw data into vector objects
    dataPoints = createVectorObjects(data)

    # create bins in grid
    grid = createGrid(nrBinsX, nrBinsY, nrBinsZ, xMin, xMax, yMin, yMax, zMin, zMax)

    # assign vector objects to correct bins
    # grid is the 3D array filled with gridBin objects containing
    # the correct vector objects
    grid = assignVectorsToGrid(dataPoints, grid)

    # report to user
    t2 = time.time()
    print("Total time: ", "{:.2f}".format(t2 - t1), " s")

    return grid


data = getRectangularGridWithVectors(10, 10, 10)
nrBinsX,nrBinsY,nrBinsZ = 10,10,10


for i in range(nrBinsX):
    for j in range(nrBinsY):
        for k in range(nrBinsZ):
            # normal averaging method
            data[i][j][k].calculateNormalAverage()

            # Gaussian averaging method
            data[i][j][k].calculateStandardDeviation()
            datasdU, datasdV, datasdW = data[i][j][k].calculateVariance()
