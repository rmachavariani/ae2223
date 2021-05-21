import gridConstructionSphereFast as gr
from class_def import * # class_def IMPORTED !!!
import numpy as np
import EnsemblePolyfit as ens

import time

pitch = [7.5]
radius = [15]

class aveObject:
    def __init__(self, inputAve):
        self.polyAve = inputAve
def saveGrid(nrOfParticles, minParticlesForAverages, pitch, radius):
    # load the grids
    grids = gr.allgrid(pitch,radius,nrOfParticles)

    t1 = time.time()

    # loop over grids to calculate averages
    for grid in grids:

        limitArray = gridBin.minmax
        dataArray = np.empty((np.size(grid, axis=0), np.size(grid, axis=1), np.size(grid, axis=2), 3), dtype=float)


        # loop over current grid

        for i in range(np.size(grid,axis=0)):
            for j in range(np.size(grid,axis=1)):
                for k in range(np.size(grid,axis=2)):

                    # current bin
                    thisBin = grid[i, j, k]

                    if len(thisBin.vectors) > minParticlesForAverages:
                        # do the polyfit
                        ens.solve(thisBin)

                        dataArray[i,j,k, :] = np.array(thisBin.polyfitAverage)
            print('i', i)
        np.save('ensaveFull', dataArray)
        np.save('limitsFull', limitArray)
        print('done saving')
#saveGrid(None, 50, pitch, radius)

############################ USING THE SAVED GRID ##################################



from math import ceil, floor, sqrt, cos, asin
import h5py

velocitydata = np.load('ensave.npy')


def determineMaxMin():
    '''Determines the mininum and maximum value for every dimension
    :return: xMin, xMax, yMin, yMax, zMin, zMax
    '''

    t1 = time.time()

    limitdata = list(np.load('limits.npy'))

    # report to user
    t2 = time.time()
    print("Max and min found in ", "{:.2f}".format(t2 - t1), " s")

    return limitdata





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
                grid[i, j, k] = gridBin(x[i], y[j], z[k], i, j, k)
                #print(velocitydata[i,j,k])
                ref = np.array([0.,0.,0.])
                comparison = ref == velocitydata[i,j,k]
                if not comparison.all():
                    grid[i,j,k].polyfitAverage = list(velocitydata[i,j,k])

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



def checkRadiusLargeEnough(pitch,radius):

    # calculate radius in between centers
    R = radius * cos(asin(pitch / (2 * radius)))

    # positive if radius is big enough
    if R >= sqrt(2)/2*pitch:
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
        #data = loadData(nrRows)

        # determine max and min of data in every dimension

        # set parameters for bins
        xMin = minMax[0]
        xMax = minMax[1]
        yMin = minMax[2]
        yMax = minMax[3]
        zMin = minMax[4]
        zMax = minMax[5]

        # transform raw data into vector objects
        #dataPoints = createVectorObjects(data)

        # create bins in grid
        grid = createGridPitchAndRadius(pitch,radius,xMin,xMax,yMin,yMax,zMin,zMax)

        # assign vector objects to correct bins
        # grid is the 3D array filled with gridBin objects containing
        # the correct vector objects
        #grid = assignVectorsToGrid(dataPoints,grid,pitch,radius,xMin,yMin,zMin)

        # report to user
        t2 = time.time()
        print("Total time: ","{:.2f}".format(t2-t1)," s")

        return grid

    else:
        print("WARNING: Set a bigger radius")


def allgrid(pitches,radii):

    t1 = time.time()

    contall = []
    for i in range(len(pitches)):
        contall.append(checkRadiusLargeEnough(pitches[i], radii[i]))

    if all(contall):

        # load in data
        minMax = determineMaxMin()

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
            #grid = assignVectorsToGrid(dataPoints, grid, pitches[i], radii[i], xMin, yMin, zMin)

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

def getGrid():
    return allgrid(pitch, radius)


#a = getGrid()
#print(a[0][25,25,25].polyfitAverage)