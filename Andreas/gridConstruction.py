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
    data = np.loadtxt("carMirrorData.dat")
    t2 = time.time()
    print("Loading done in ","{:.2f}".format(t2-t1)," s")
    
    return data

def determineMaxMin(data):
    '''Determines the mininum and maximum value for every dimension
    :return: xMin, xMax, yMin, yMax, zMin, zMax
    '''

    t1 = time.time()

    # determine min and max
    xMin = np.amin(data[:,0])
    xMax = np.amax(data[:,0])
    yMin = np.amin(data[:,1])
    yMax = np.amax(data[:,1])
    zMin = np.amin(data[:,2])
    zMax = np.amax(data[:,2])

    # report to user
    t2 = time.time()
    print("Max and min found in ","{:.2f}".format(t2-t1)," s")

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
    print("Objects created in ","{:.2f}".format(t2-t1)," s")

    return dataPoints

def createGrid(nrBinsX,nrBinsY,nrBinsZ,xMin,xMax,yMin,yMax,zMin,zMax):
    '''Creates the grid by generating a 3D numpy array filled with
    objects of class gridbin
    :return: nrX x nrY x nrZ numpy array
    '''

    t1 = time.time()
    
    # create empty 3D array
    grid = np.empty((nrBinsX,nrBinsY,nrBinsZ),dtype=object)
    
    # calculate width of bins in all directions
    widthX = (xMax-xMin) / nrBinsX
    widthY = (yMax-yMin) / nrBinsY
    widthZ = (zMax-zMin) / nrBinsZ

    # set widths of bin to class static members
    gridBin.widthX = widthX
    gridBin.widthY = widthX
    gridBin.widthZ = widthZ
    
    # define x, y and z coordinates of center bin
    x = np.linspace(xMin,xMax-widthX,nrBinsX) + widthX/2
    y = np.linspace(yMin,yMax-widthY,nrBinsY) + widthY/2
    z = np.linspace(zMin,zMax-widthZ,nrBinsZ) + widthZ/2
    
    # fill matrix with bin objects by looping over matrix
    for i in range(nrBinsX):
        for j in range(nrBinsY):
            for k in range(nrBinsZ):
                
                grid[i,j,k] = gridBin(x[i],y[j],z[k])

   # report to user
    t2 = time.time()
    print('Grid created in ',"{:.2f}".format(t2-t1)," s")
   
    return grid


def assignVectorsToGrid(vectors,grid):
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
        iX = floor( x / widthX )
        iY = floor( y / widthY )
        iZ = floor( z / widthZ )
        
        # assign to correct bin
        grid[iX,iY,iZ].addVector(vector)
        
    # report to user
    t2 = time.time()
    print("Assigning of vectors to bins completed in ","{:.2f}".format(t2-t1)," s")
        
    return grid

def showAmountOfVectorsInBin(grid):
    ''' reports on amount of vectors in every bin'''

    sizeX = np.size(grid,axis=0)
    sizeY = np.size(grid,axis=1)
    sizeZ = np.size(grid,axis=2)
    
    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                
                print(len(grid[i][j][k].vectors))





#------------------------------MAIN-----------------------------#

def getGridWithVectors(nrBinsX,nrBinsY,nrBinsZ):

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
    grid = createGrid(nrBinsX,nrBinsY,nrBinsZ,xMin,xMax,yMin,yMax,zMin,zMax)

    # assign vector objects to correct bins
    # grid is the 3D array filled with gridBin objects containing
    # the correct vector objects
    grid = assignVectorsToGrid(dataPoints,grid)

    # report to user
    t2 = time.time()
    print("Total time: ","{:.2f}".format(t2-t1)," s")

    return grid

data = getGridWithVectors(10,10,10)


