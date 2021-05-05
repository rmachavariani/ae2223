import numpy as np
from class_def import *
from math import ceil, floor, sqrt, cos, asin
import time
import h5py


def loadData(nrRows):
    ''' Loads data from carMirrorData.dat
    :return: nrVectors x 6 np array
    '''
    # load the data
    t1 = time.time()

    #use for 2.7 million particles
    #data = np.loadtxt("/Users/stefanrooze/Downloads/dataset.txt",max_rows = nrRows)

    #use for 91 million particles
    raw_data = h5py.File('/Users/stefanrooze/Downloads/completeDatasetCarMirrorTracksNikhilesh.dat', 'r')
    data_column = raw_data['output']
    data_column = np.array(data_column)
    data = np.transpose(data_column)
    print(data[100])
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


def checkRadiusLargeEnough(pitch,radius):

    # calculate radius in between centers
    #R = radius * cos(asin(pitch / (2 * radius)))
    dis = sqrt(pitch**2)
    # positive if radius is big enough
    if dis <= 2*radius:#R >= sqrt(2)/2*pitch:
        return True
    else:
        return False



#-------------------------------MAIN--------------------------------#




def getSphericalGridWithVectorsFast(pitch,radius,nrRows):
    '''

    :param pitch: pitch between centers of bins in mm
    :param radius: radius of bins in mm
    :param nrRows: amount of rows to load from the datafile
                    enter None when all have to be loaded
    :return: grid with bins and assigned objects
    '''

    t1 = time.time()

    # check if radius is big enough
    cont = checkRadiusLargeEnough(pitch,radius)

    # only perform creation of grid when every particle is covered
    if cont:

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
        grid = createGridPitchAndRadius(pitch,radius,xMin,xMax,yMin,yMax,zMin,zMax)

        # assign vector objects to correct bins
        # grid is the 3D array filled with gridBin objects containing
        # the correct vector objects
        grid = assignVectorsToGrid(dataPoints,grid,pitch,radius,xMin,yMin,zMin)

        # report to user
        t2 = time.time()
        print("Total time: ","{:.2f}".format(t2-t1)," s")

        return grid

    else:
        print("WARNING: Set a bigger radius")


def loadParticles(nrRows):
    # load data
    data = loadData(nrRows)

    # transform raw data into vector objects
    dataPoints = createVectorObjects(data)

    # determine min and max
    minMax = determineMaxMin(data)

    return dataPoints, minMax

def allgrid(pitches,radii,nrParticles):

    t1 = time.time()

    contall = []
    for i in range(len(pitches)):
        contall.append(checkRadiusLargeEnough(pitches[i], radii[i]))

    if all(contall):

        # load in data
        dataPoints, minMax = loadParticles(nrParticles)

        # set parameters for bins
        xMin = minMax[0]
        xMax = minMax[1]
        yMin = minMax[2]
        yMax = minMax[3]
        zMin = minMax[4]
        zMax = minMax[5]

        grids = []
        # create different grids
        for i in range(len(pitches)):

            print()
            print("---------------")
            print("Grid ", i+1,)


            # create the grid with certain pitch and radius
            grid = createGridPitchAndRadius(pitches[i], radii[i], xMin, xMax, yMin, yMax, zMin, zMax)

            # assign vector objects to correct bins
            # grid is the 3D array filled with gridBin objects containing
            # the correct vector objects
            grid = assignVectorsToGrid(dataPoints, grid, pitches[i], radii[i], xMin, yMin, zMin)

            grids.append(grid)

        # report to user
        t2 = time.time()
        print("---------------")
        print()
        print("Total time: ", "{:.2f}".format(t2 - t1), " s")
        print()
        print("----------------------------------------")

        return grids

    else:
            print("WARNING: Set a bigger radius")