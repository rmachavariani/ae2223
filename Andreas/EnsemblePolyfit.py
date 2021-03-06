import numpy as np
import pdb
import matplotlib.pyplot as plt
import matplotlib
from functools import partial as part
import statistics as sts
matplotlib.use("Qt5Agg")


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

# Calculates angular velocity of the vorticity field


def vorticity(ucoefs, vcoefs, wcoefs):
    curl = np.empty(3)
    curl[0] = (wcoefs[2] - vcoefs[3]) / 2
    curl[1] = (ucoefs[3] - wcoefs[1]) / 2
    curl[2] = (vcoefs[1] - ucoefs[2]) / 2

    return curl


def solve(bin):
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
    bin.polyfitAverage = [ucoefs[0], vcoefs[0], wcoefs[0]]
    bin.vorticity = vorticity(ucoefs, vcoefs, wcoefs)


# data = getRectangularGridWithVectors(10, 10, 10)
# test_cell = data[5, 5, 5]
# solve(test_cell)
# print(test_cell.polyfitAverage)
# print(test_cell.vorticity)
# print(test_cell.fluc)
# print(test_cell.turb_eng)



