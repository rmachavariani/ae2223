import class_def
import numpy as np
from gridConstructionSphereFast import *
from gridConstruction import *
def get_nodes_array(location, scheme, direction):
    scheme_dict_shift = {'forward': 0, 'center': -2, 'backward': -4}

    shift = np.zeros((5, 3))
    for i in range(np.size(shift, axis=0)):
        for j in range(np.size(shift, axis=1)):
            if j == direction:
                shift[i, j] = i + scheme_dict_shift[scheme]

    for i in range(np.size(shift, axis=0)):
        shift[i,:] += location

    return shift.astype(int)

def partial_approx(grid, scheme, veldir, griddir,poss,h):

    scheme_weights = {'forward': np.array([-25,48,-36,16,-3]), 'central' : np.array([1,-8,0,8,-1]), 'backward' : np.array([3,-16,36,-48,25])}

    nodes = get_nodes_array(poss, scheme, griddir)

    weighted_node_sum = 0

    for i in range(np.size(nodes, axis=0)):
        print(nodes[i,1])
        #print(grid[nodes[i,1]])
        print(range(np.size(nodes,axis=0)))
        i = nodes[i,0]
        j = nodes[i,1]
        k = nodes[i,2]
        weighted_node_sum += grid[i,j,k].polyfitAverage[veldir]*scheme_weights[scheme][i]

    return weighted_node_sum / (12 * h)


def vorticity(grid, bin):
    # determine step size [m]
    h = (grid[0, 0, 0].x - grid[1, 0, 0].x) / 1000

    # Initialize the curl definition
    pdes = ['dwdy', 'dvdz', 'dudz', 'dwdx', 'dvdx', 'dudy']
    veldir = [2, 1, 0, 2, 1, 0]     # velocity components in vorticity defintion ()
    griddir = [0, 2, 2, 0, 0, 1]

    # Get position of the bin
    poss = np.array([bin.i, bin.j, bin.k])
    limits = [gridBin.nrBinsX - 1, gridBin.nrBinsY - 1, gridBin.nrBinsZ - 1]
    idx = 0
    scheme = None

    # Determine scheme for the partial derivative
    for pd in griddir:

        if poss[pd] < 2:
            scheme = 'forward'
        elif poss[pd] > (limits[pd] - 2):
            scheme = 'backward'
        else:
            scheme = 'central'

        pdes[idx] = partial_approx(grid, scheme, veldir, griddir, poss, h)
    idx += 1

    return [pdes[0] - pdes[1], pdes[2] - pdes[3], pdes[4] - pdes[5]]


# DEBUGGING:
testGrid = getRectangularGridWithVectors(15,15,15)
print(testGrid[5,5,5])
print(vorticity(testGrid,testGrid[5,5,5]))

