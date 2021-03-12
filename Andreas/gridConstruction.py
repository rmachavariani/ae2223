import numpy as np
from class_def import *
from math import floor
import time


def loadData():
    
    # load the data 
    t1 = time.time()
    data = np.loadtxt("carMirrorData.dat")
    t2 = time.time()
    print("Loading done")
    print(t2-t1)
    
    return data
    
    
def createGrid(nrBinsX,nrBinsY,nrBinsZ,xMin,xMax,yMin,yMax,zMin,zMax):
    
    t1 = time.time()
    
    # create empty 3D array
    grid = np.empty((nrBinsX,nrBinsY,nrBinsZ),dtype=object)
    
    # calculate width of bins in all directions
    widthX = (xMax-xMin) / nrBinsX
    widthY = (yMax-yMin) / nrBinsY
    widthZ = (zMax-zMin) / nrBinsZ
    
    gridBin.widthX = widthX
    gridBin.widthY = widthX
    gridBin.widthZ = widthZ
    
    # define x, y and z points
    x = np.linspace(xMin,xMax-widthX,nrBinsX)
    y = np.linspace(yMin,yMax-widthY,nrBinsY)
    z = np.linspace(zMin,zMax-widthZ,nrBinsZ)
    
    # fill matrix with bin objects
    for i in range(nrBinsX):
        for j in range(nrBinsY):
            for k in range(nrBinsZ):
                
                grid[i,j,k] = gridBin(x[i],y[j],z[k])
   
    t2 = time.time()
    print('Grid created in ', t2-t1)
   
    return grid


def assignVectorsToGrid(vectors,grid):
    
    # start time
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
        
    # stop time and print message
    t2 = time.time()
    print("Assigning of vectors to bins completed")
    print("Total time: ",t2-t1," s")
        
    return grid

def showAmountOfVectorsInBin(grid):
    
    sizeX = np.size(grid,axis=0)
    sizeY = np.size(grid,axis=1)
    sizeZ = np.size(grid,axis=2)
    
    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                
                print(len(grid[i][j][k].vectors))

def determineMaxMin(data):

    t1 = time.time()

    xMin = np.amin(data[:,0])
    xMax = np.amax(data[:,0])
    yMin = np.amin(data[:,1])
    yMax = np.amax(data[:,1])
    zMin = np.amin(data[:,2])
    zMax = np.amax(data[:,2])

    t2 = time.time()
    print("Max and min found")
    print("Duration: ",t2-t1)

    return xMin, xMax, yMin, yMax, zMin, zMax

def createVectorObjects(data):
    
    t1 = time.time()
    dataPoints = np.empty(np.size(data,axis=0),dtype=object)
    for i in range(np.size(data,axis=0)):
        dataPoints[i] = vector(data[i,:])
        
    t2 = time.time()
    print("Objects created")
    print(t2-t1)
    
    return dataPoints

#------------------------------MAIN-----------------------------#

t1 = time.time()

data = loadData()
minMax = determineMaxMin(data)

nrBinsX = 10
nrBinsY = 10
nrBinsZ = 10
xMin = minMax[0]
xMax = minMax[1]
yMin = minMax[2]
yMax = minMax[3]
zMin = minMax[4]
zMax = minMax[5]


dataPoints = createVectorObjects(data)

grid = createGrid(nrBinsX,nrBinsY,nrBinsZ,xMin,xMax,yMin,yMax,zMin,zMax)

grid = assignVectorsToGrid(dataPoints,grid)

showAmountOfVectorsInBin(grid)

t2 = time.time()
print("Total time: ", t2-t1)




