import gridConstructionSphereFast as gr
from class_def import * # class_def IMPORTED !!!
import numpy as np
import EnsemblePolyfit as ens
import vorticity as vort
import time
import ArrayGenerator as ag

# parameters --> SET THESE IN ARRAYGENERATOR
# nrOfParticles = None
# minParticlesForAverages = 50
#
# pitch = [10]
# radius = [10]

# load the grids
grids = ag.getGrid()


t1 = time.time()

# loop over grids to calculate averages





print('Commencing vorticity calculation')
for grid in grids:

    # loop over current grid
    for i in range(np.size(grid,axis=0)):
        for j in range(np.size(grid,axis=1)):
            for k in range(np.size(grid,axis=2)):

                # current bin
                thisBin = grid[i, j, k]

                if len(thisBin.vectors) > minParticlesForAverages:
                    thisBin.vorticity = vort.vorticity(grid, thisBin, [gridBin.nrBinsX - 1, gridBin.nrBinsY - 1, gridBin.nrBinsZ - 1], 3)
                    #print(thisBin.vorticity)

# TESTING::


#print(grids[0][50,30,25].vorticity)