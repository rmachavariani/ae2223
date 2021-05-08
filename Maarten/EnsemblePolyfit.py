#from gridConstruction import *
import pdb
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from functools import partial as part
import statistics as sts
matplotlib.use("Qt5Agg")
import Vorticity as vort

def fetch_vector(bin):
    test_vectors = bin.vectors
    particles_val = np.empty((len(test_vectors), 6))

    if test_vectors == [] or () or None:
        print("No particle")
    else:
        for i in range(len(test_vectors)):
            x = test_vectors[i].x
            y = test_vectors[i].y
            z = test_vectors[i].z
            u = test_vectors[i].u
            v = test_vectors[i].v
            w = test_vectors[i].w

            particles_val[i] = [x, y, z, u, v, w]
    # print(particles_val[:, 3])
    return particles_val


def basis(particles_val, cell):
    design_matrix = np.empty((len(particles_val), 10))
    for i in range(len(particles_val)):
        dx = (particles_val[i, 0] - cell.x) / 1000
        dy = (particles_val[i, 1] - cell.y) / 1000
        dz = (particles_val[i, 2] - cell.z) / 1000

        design_matrix[i] = [1, dx, dy, dz, dx * dy, dx * dz, dy * dz, dx ** 2, dy ** 2, dz ** 2]

    # print(design_matrix)
    return design_matrix


def coefficients(basis, fnc):
    coeffs = np.linalg.lstsq(basis, fnc, rcond=None)[0]

    return coeffs



def createPolyFit(coefs):

    # create partical function
    partialFit = part(polyFit, coefs)

    return partialFit


def polyFit(coefs, dx, dy, dz):

    # define basis
    basis = np.array([1, dx, dy, dz, dx * dy, dx * dz, dy * dz, dx ** 2, dy ** 2, dz ** 2])

    # calculate value
    functionValue = np.sum(basis * coefs)

    return functionValue


def vorticity(ucoefs, vcoefs, wcoefs):
    curl = np.empty(3)
    curl[0] = wcoefs[2] - vcoefs[3]
    curl[1] = ucoefs[3] - wcoefs[1]
    curl[2] = vcoefs[1] - ucoefs[2]

    return curl






def solve(bin):
    #print('SOLVE EXECUTED')
    data = fetch_vector(bin)
    design_matrix = basis(data, bin)
    ucoefs = coefficients(design_matrix, data[:, 3])
    vcoefs = coefficients(design_matrix, data[:, 4])
    wcoefs = coefficients(design_matrix, data[:, 5])
    inst_vel = data[:, 3:6]

    avg_vel = [ucoefs[0], vcoefs[0], wcoefs[0]]
    vel_prime = [sts.stdev(inst_vel[:, 0], avg_vel[0]), sts.stdev(inst_vel[:, 1], avg_vel[1]), sts.stdev(inst_vel[:, 2], avg_vel[2])]
    energy = 0.5 * (vel_prime[0]**2 + vel_prime[1]**2 + vel_prime[2]**2)

    bin.fluc = vel_prime
    bin.turb_eng = energy
    bin.fitU = createPolyFit(ucoefs)
    bin.fitV = createPolyFit(vcoefs)
    bin.fitW = createPolyFit(wcoefs)
    #print(ucoefs[0], vcoefs[0], wcoefs[0])
    bin.polyfitAverage = [ucoefs[0], vcoefs[0], wcoefs[0]]
    #print('POLYFIT EXECUTED')
    #bin.vorticity = vort.vorticity(grid, bin, limits)


# data = getRectangularGridWithVectors(10, 10, 10)
#
# sizeX = np.size(data, axis=0)
# sizeY = np.size(data, axis=1)
# sizeZ = np.size(data, axis=2)

# for i in range(sizeX):
#     for j in range(sizeY):
#         for k in range(sizeZ):
#             cell = data[i,j,k]
#             if len(cell.vectors) > 50:
#                 solve(cell)
#                 print('---------')
#                 counter  = (str(i) + ',' + str(j) + ',' + str(k))
#                 print('particles in bin: ' + str(len(cell.vectors)))
#                 print('cell #: ' + counter)
#                 print('particles in bin: ')
#                 print('cell x: ' + str(cell.x) + ' cell y: ' + str(cell.y) + ' cell z: '+ str(cell.z))
#                 print('u: ' + str(cell.polyfitAverage[0]) + ' v: ' + str(cell.polyfitAverage[1]) + ' w: ' +str(cell.polyfitAverage[2]))
#                 u, v, w = cell.polyfitAverage[0], cell.polyfitAverage[1], cell.polyfitAverage[2]
#                 if u > 15 or u < -5 or abs(v) > 5 or abs(w) > 5:
#                     print(u, v, w)
#                     print('outlier detected!')
#
#
#
#
#                 print('---------')



#print(test_cell.polyfitAverage)
#print(test_cell.vorticity)
#print(test_cell.fluc)
#print(test_cell.turb_eng)



