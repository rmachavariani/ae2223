import gridConstructionSphereFast as gr
import numpy as np
import EnsemblePolyfit as ens
import time

# parameters
nrOfParticles = None

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
                thisBin = grid[i,j,k]

                # calculate normal average
                #thisBin.calculateNormalAverage()

                # calculate variance and deviation
                #thisBin.calculateStandardDeviation()
                #thisBin.calculateVariance()

                # calculate gaussian average
                #thisBin.calculateGaussianAverage()

                # calculate dog gaussian average


                # calculate poly fit
                if len(thisBin.vectors) > 10:
                    ens.solve(thisBin)

t2 = time.time()

print()
print("Total time for taking averages: ", round(t2-t1,2))
print()
print("----------------------------------------")

