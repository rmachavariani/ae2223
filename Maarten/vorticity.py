import class_def
from gridConstructionSphereFast import *

def vorticity(grid, bin):

    # lists to loop over the bin
    bin_loc = [bin.i, bin.j, bin.k]
    grid_size = [gridBin.nrBinsX, gridBin.nrBinsY, gridBin.nrBinsZ]
    pdus = [[],[],[]]

    h = grid[0,0,0].x - grid[1,0,0].x

    for cmp in bin_loc:

        if cmp < 2 or cmp > (grid_size[bin_loc.index(cmp)]-2):
            partial_derivative1 = edge_scheme()


        else:
            central_scheme()



def calc_vorticity(grid, bin):
    # initialize the six pd's
    h = grid[0, 0, 0].x - grid[1, 0, 0].x
    pdes = [dwdy, dvdz, dudz, dwdx, dvdx, dudy]
    veldir = [2,1,0,2,1,0]
    griddir = [0, 2, 2, 0, 0, 1]
    poss = np.array([bin.i, bin.j, bin.k])
    limits = [gridBin.nrBinsX - 1, gridBin.nrBinsY - 1, gridBin.nrBinsZ -1]
    idx = 0
    for pd in griddir:

        if poss[pd] < 2:



            shift = np.zeros(4,3)
            for i in range(np.size(shift, axis=0)):
                for j in range(np.size(shift, axis=1)):
                    if j == pd:
                        shift[i,j] = i+1
            umin = bin.polyfitAverage[veldir[idx]]
            u1 =
            jcor =
            kcor =
            pdes[idx] = edge_scheme(grid, pd, veldir[idx] poss, h)

        elif poss[pd] > (limits[pd] - 2):


        else:
            central_scheme()

        idx += 1




    # calculate pd1 for vcomp --> identify suitable scheme
    # calculate pd2 for vcomp --> identify suitable scheme
    # return the vorticity

    # direction of partial derivative: griddir
    # grid[i,j,k]
    # griddir[idx]: varying
    # other values fixed
    # np.zeros(5,3)


def central_scheme(grid, direction, component, location, h):
    nodes = [node_1, node_2, node_3, node_4, node_5]
    shift = np.zeros(5,3)
    for i in range(np.size(shift, axis=0)):
        for j in range(np.size(shift, axis=1)):
            if j == direction:
                shift[i,j] = j

    for rows in shift:
       nodes[rows] = shift[i,:] + position
    umin2 = location[0] + shift[0,0]
    umin1 =
    u1 =
    u2 =

    return (1*umin2-8*umin1+8*u1-1*u2)/(12.0*h)


def left_scheme(grid, direction, location, h):
    uzero =
    u1 =
    u2 =
    u3 =
    u4 =
    return (-25*uzero+48*u1-36*u2+16*u3-3*u4)/(12.0*h)

def right_scheme(grid, direction, location, h):
    umin4 =
    umin3 =
    umin2 =
    umin1 =
    uzero =

    return (3*umin4-16*umin3+36*umin2-48*umin1+25*uzero)/(12.0*h)