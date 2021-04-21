import gridConstructionSphereFast as gr
import numpy as np
import EnsemblePolyfit as ens
import time

# parameters
nrOfParticles = None
minParticlesForAverages = 50

pitch = [10]
radius = [10]

# load the grids
grids = gr.allgrid(pitch,radius,nrOfParticles)


t1 = time.time()

# loop over grids to calculate averages
for grid in grids:

    # loop over current grid
    for i in range(np.size(grid,axis=0)):
        for j in range(np.size(grid,axis=1)):
            for k in range(np.size(grid,axis=2)):

                # current bin
                thisBin = grid[i, j, k]

                if len(thisBin.vectors) > minParticlesForAverages:

                    # calculate normal/top hat average
                    thisBin.calculateNormalAverage()

                    # calculate variance and deviation
                    # for normal/top hat average
                    thisBin.calculateStandardDeviation()
                    thisBin.calculateVariance()

                    # calculate gaussian average
                    thisBin.calculateGaussianAverage()

                    # calculate poly fit
                    ens.solve(thisBin)

t2 = time.time()

print()
print("Total time for taking averages: ", round(t2-t1,2))
print()
print("----------------------------------------")


