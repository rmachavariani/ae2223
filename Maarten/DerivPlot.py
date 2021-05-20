import matplotlib.pyplot as plt
import vorticity as vort
import gridConstructionSphereFast as gr
import time
import numpy as np
import EnsemblePolyfit as ens
import class_def
import pdb
import ArrayGenerator as ag


# parameters
nrOfParticles = None
minParticlesForAverages = 50


pitch = [7.5]
radius = [15]

# load the grids
grids = ag.getGrid()


t1 = time.time()

# loop over grids to calculate averages






#############################################################

input_z = 75
input_x = 75
# Find the i corresponding to x and the k corresponding to z
x = 0
z = 0

for i in range(np.size(grids[0], axis=0)):
    if abs(grids[0][i, 0, 0].x - input_x) < radius[0]:
        x = i

# Find index for bin with bin.z closest to z specified --> k

for k in range(np.size(grids[0], axis=2)):
    if abs(grids[0][0, 0, k].z - input_z) < radius[0]:
        z = k

# get plot limits: ymin, ymax
# ymin = grids[0][0, 0, 0].y
# ymax = grids[0][0, z, 0]
# Make sure the partial derivatives are stored somewhere DONE

print('Commencing vorticity calculation')
for grid in grids:
    # loop over current grid
    for j in range(np.size(grid,axis=1)):
        thisBin = grid[x, j, z]
        #print(thisBin.polyfitAverage)
        if thisBin.polyfitAverage:
            thisBin.vorticity = vort.vorticity(grid, thisBin, [np.size(grid, axis=0) - 1, np.size(grid, axis=1) - 1, np.size(grid, axis=2) - 1], 3)



deriv_uy = []
deriv_uz = []
y_coords = []
z_coords = []

for j in range(np.size(grid, axis=1)):
    grid = grids[0]
    cbin = grid[x,j,z]
    if cbin.pde:
        pd = cbin.pde[0]
        y_coords.append(cbin.y)
        deriv_uy.append(pd)


    #deriv_uz.append(grids[0][x, j, z].pde[1])
    #for elements
    #y_coords.append(grids[0][x, j, z].y)
    #z_coords.append(grids[0][x, j, z].z)
#print(deriv_uy)
#print(y_coords)

plt.plot(y_coords, deriv_uy)
plt.xlabel('y position', fontsize=15)
plt.ylabel("du/dy", fontsize=15)
plt.grid()
plt.show()

# do same for x --> i
# for j in range(the length of j)
# plot the dudy along [i,j,k]
# plot(dudy(grid(i,j,k).pde[0]), grid(i,j,k).y)

# create array with