import numpy as np
from class_def import *
from math import ceil, floor, sqrt, cos, asin
import time
import matplotlib.pyplot as plt


def loadData(nrRows):
    ''' Loads data from carMirrorData.dat
    :return: nrVectors x 6 np array
    '''
    # load the data
    t1 = time.time()
    data = np.loadtxt("/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/My coding/carMirrorData.dat",max_rows = nrRows)
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

    # report amount of bins to user
    xAmount = np.size(grid, axis=0)
    yAmount = np.size(grid, axis=1)
    zAmount = np.size(grid, axis=2)
    print("Amount of bins in x direction: ", xAmount)
    print("Amount of bins in y direction: ", yAmount)
    print("Amount of bins in z direction: ", zAmount)
    print("Total amount of bins: ", xAmount * yAmount * zAmount)

    return grid,xAmount,yAmount,zAmount

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

    # report amount of bins to user
    xAmount = np.size(grid,axis=0)
    yAmount = np.size(grid,axis=1)
    zAmount = np.size(grid,axis=2)
    print("Amount of bins in x direction: ", xAmount)
    print("Amount of bins in y direction: ",yAmount)
    print("Amount of bins in z direction: ",zAmount)
    print("Total amount of bins: ",xAmount*yAmount*zAmount)

    return grid


def assignVectorsToGrid(vectors, grid, pitch, radius,
                        xMin, yMin, zMin):

    t1 = time.time()

    # loop though all the vectors
    for vector in vectors:

        # get coordinates
        x = vector.x
        y = vector.y
        z = vector.z

        # calculate indices in every direction
        indexXLow = int(ceil((x-radius-xMin) / pitch))
        indexXHigh = int(floor((x + radius - xMin) / pitch))
        indexYLow = int(ceil((y - radius - yMin) / pitch))
        indexYHigh = int(floor((y + radius - yMin) / pitch))
        indexZLow = int(ceil((z - radius - zMin) / pitch))
        indexZHigh = int(floor((z + radius - zMin) / pitch))

        # create range of indices in every direction
        xRange = range(indexXLow,indexXHigh + 1)
        yRange = range(indexYLow,indexYHigh + 1)
        zRange = range(indexZLow,indexZHigh + 1)

        # loop through all relevant bins
        for i in xRange:
            for j in yRange:
                for k in zRange:

                    # get bin object
                    aBin = grid[i][j][k]

                    # get coordinates
                    xx = aBin.x
                    yy = aBin.y
                    zz = aBin.z

                    if sqrt((xx-x)**2 + (yy-y)**2 + (zz-z)**2) <= radius:

                        aBin.addVector(vector)

    # report to user
    t2 = time.time()
    print("Assigning of vectors to bins completed in ", "{:.2f}".format(t2 - t1), " s")

    return grid

def showBins(pitch,radius):

    edge = 5
    outermost = 2 * radius + pitch + 2 * edge

    try:
        # calculate new radius
        radius = radius * cos(asin(pitch / (2 * radius)))


        outermost = 2 * radius + pitch + 2 * edge
        circle1 = plt.Circle((radius + edge, radius + edge),radius,edgecolor="black")
        circle2 = plt.Circle((radius + pitch + edge,radius + edge),radius,edgecolor="black")
        circle3 = plt.Circle((radius + edge,radius + pitch + edge),radius,edgecolor="black")
        circle4 = plt.Circle((radius + pitch + edge,radius + pitch + edge),radius,edgecolor="black")

        fig, ax = plt.subplots()  # note we must use plt.subplots, not plt.subplot
        # (or if you have an existing figure)
        # fig = plt.gcf()
        # ax = fig.gca()

        ax.add_patch(circle1)
        ax.add_patch(circle2)
        ax.add_patch(circle3)
        ax.add_patch(circle4)

        ax.set_xlim([0,outermost])
        ax.set_ylim([0,outermost])

        ax.set_xticks([])
        ax.set_yticks([])


        plt.show()

    except:
        print("Test")
        fig, ax = plt.subplots()
        ax.set_xlim([0,outermost])
        ax.set_ylim([0,outermost])

        ax.set_xticks([])
        ax.set_yticks([])


        plt.show()


#-------------------------------MAIN--------------------------------#

def getSphericalGridWithVectors(pitch,radius,nrRows):

    t1 = time.time()

    # show grid
    showBins(pitch,radius)
    user = input("Continue? (y/n): ")

    if user == "y":
        # load the data
        data = loadData(nrRows)

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
        grid,xAmount,yAmount,zAmount = createGridPitchAndRadius(pitch,radius,xMin,xMax,yMin,yMax,zMin,zMax)

        # assign vector objects to correct bins
        # grid is the 3D array filled with gridBin objects containing
        # the correct vector objects
        grid = assignVectorsToGrid(dataPoints,grid,pitch,radius,xMin,yMin,zMin)

        # report to user
        t2 = time.time()
        print("Total time: ","{:.2f}".format(t2-t1)," s")

        return grid,xAmount,yAmount,zAmount

    else:

        print("Grid creation was not performed")

grid,xAmount,yAmount,zAmount = getSphericalGridWithVectors(30,30,1000)

for i in range(xAmount):
    for j in range(yAmount):
        for k in range(zAmount):
            # normal averaging method
            grid[i][j][k].calculateNormalAverage()

            # Gaussian averaging method
            grid[i][j][k].calculateStandardDeviation()
            grid[i][j][k].calculateGaussianAverage()
