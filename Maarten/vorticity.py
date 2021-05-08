#   Using vorticity:
#   1. See Main.py in my folder, an extra loop is made: this is essential because
#   all bin averaging has to be completed before the vorticity can be calculated
#
#   2. the class definition has the i, j and k added to gridBin, this must also be added in the
#   gridConstructionSphereFast function were the instances of gridbin are placed in the grid
#
#   3. class_def imported in Main.py -> needed for the argument of vorticity()

# ========================================================
import class_def
import numpy as np
import gridConstructionSphereFast as gr
#from gridConstruction import *
import EnsemblePolyfit as ens
def get_nodes_array(location, scheme, direction):
    #print(location, scheme, direction)
    scheme_dict_shift = {'forward': 0, 'central': -2, 'backward': -4}

    shift = np.zeros((5, 3))
    for i in range(np.size(shift, axis=0)):
        for j in range(np.size(shift, axis=1)):
            if j == direction:
                shift[i, j] = i + scheme_dict_shift[scheme]
    #print(shift)
    for i in range(np.size(shift, axis=0)):
        shift[i,:] += location
    #print(shift)
    return shift.astype(int)

def partial_approx(grid, scheme, veldir, griddir,poss,h, noData):

    scheme_weights = {'forward': np.array([-25,48,-36,16,-3]), 'central' : np.array([1,-8,0,8,-1]), 'backward' : np.array([3,-16,36,-48,25])}

    nodes = get_nodes_array(poss, scheme, griddir)
    #print(nodes)
    weighted_node_sum = 0
    for cntr in range(np.size(nodes, axis=0)):

        i = nodes[cntr,0]
        j = nodes[cntr,1]
        k = nodes[cntr,2]


        node_bin = grid[i,j,k]
        if not node_bin.polyfitAverage:
            noData = True

        else:
            node_val = node_bin.polyfitAverage[veldir]
            node_weight = scheme_weights[scheme][cntr]
            weighted_node_sum += node_val * node_weight

        if not noData:
            return weighted_node_sum / (12 * h), noData
        else:
            return 1, noData


def vorticity(grid, bin, limits):
    # determine step size [m]
    h = (grid[0, 0, 0].x - grid[1, 0, 0].x) / 1000

    # Initialize the curl definition
    pdes = ['dwdy', 'dvdz', 'dudz', 'dwdx', 'dvdx', 'dudy']
    veldir = [2, 1, 0, 2, 1, 0]     # velocity components in vorticity defintion ()
    griddir = [0, 2, 2, 0, 0, 1]

    # Get position of the bin
    poss = np.array([bin.i, bin.j, bin.k])

    idx = 0
    scheme = None
    noData = False

    # Determine scheme for the partial derivative
    for pd in griddir:

        if poss[pd] < 2:
            scheme = 'forward'
        elif poss[pd] > (limits[pd] - 2):
            scheme = 'backward'
        else:
            scheme = 'central'

        pdes[idx], noData = partial_approx(grid, scheme, veldir[idx], pd, poss, h, noData)
        idx += 1
    if not noData:
        return [pdes[0] - pdes[1], pdes[2] - pdes[3], pdes[4] - pdes[5]]
    else:
        return []


# DEBUGGING:
#testGrid = getRectangularGridWithVectors(5,5,5)

#print(testGrid[5,5,5])
#print(vorticity(testGrid,testGrid[3,2,3]))

#test = get_nodes_array([5,5,5], 'central', 1)
#print(test)

# ===============================================================

# parameters
# nrOfParticles = None
#
# pitch = [18.081]    # [10,15,20]
# radius = [13]   # [10,15,20]
#
# # load the grids
# grids = gr.allgrid(pitch, radius, nrOfParticles)
# print(grids)
#
# t1 = time.time()
#
# # loop over grids to calculate averages
# for grid in grids:
#
#     # loop over current grid
#     for i in range(np.size(grid, axis=0)):
#         for j in range(np.size(grid, axis=1)):
#             for k in range(np.size(grid, axis=2)):
#
#                 # current bin
#                 thisBin = grid[i, j, k]
#
#                 # calculate normal average
#                 # thisBin.calculateNormalAverage()
#
#                 # calculate variance and deviation
#                 # thisBin.calculateStandardDeviation()
#                 # thisBin.calculateVariance()
#
#                 # calculate gaussian average
#                 # thisBin.calculateGaussianAverage()
#
#                 # calculate dog gaussian average
#
#                 # calculate poly fit
#                 if len(thisBin.vectors) > 50:
#                     ens.solve(thisBin)
#                     thisBin.vorticity = vorticity(grid,thisBin)
#
# t2 = time.time()
#
# print()
# print("Total time for taking averages: ", round(t2-t1, 2))
# print()
# print("----------------------------------------")
