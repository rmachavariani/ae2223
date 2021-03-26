from gridConstruction import *
import pdb
import matplotlib.pyplot as plt
import matplotlib
from functools import partial as part

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
        dx = (particles_val[i, 0] - cell.x)
        dy = (particles_val[i, 1] - cell.y)
        dz = (particles_val[i, 2] - cell.z)

        design_matrix[i] = [1, dx, dy, dz, dx * dy, dx * dz, dy * dz, dx ** 2, dy ** 2, dz ** 2]

    # print(design_matrix)
    return design_matrix


def  coefficeints(basis, fnc):
    coeffs = np.linalg.lstsq(basis, fnc)[0]

    return coeffs

def solve(bin):
    data = fetch_vector(bin)
    design_matrix = basis(data, bin)
    ucoefs = coefficeints(design_matrix, data[:,3])
    vcoefs = coefficeints(design_matrix, data[:,4])
    wcoefs = coefficeints(design_matrix, data[:,5])
    bin.fitU = createPolyFit(ucoefs)
    bin.fitV = createPolyFit(vcoefs)
    bin.fitW = createPolyFit(wcoefs)
    bin.polyfitAverage.append([ucoefs[0], vcoefs[0], wcoefs[0]])
    bin.vorticity = curl(ucoefs, vcoefs, wcoefs)

def createPolyFit(coefficients):

    # create partical function
    partialFit = part(polyFit,coefficients)

    return partialFit

def polyFit(coefficients,dx,dy,dz):

    # define basis
    basis = np.array([1, dx, dy, dz, dx * dy, dx * dz, dy * dz, dx ** 2, dy ** 2, dz ** 2])

    # calculate value
    functionValue = np.sum(basis * coefficients)

    return functionValue

def curl(ucoefs, vcoefs, wcoefs):
    curl = []
    curl.append(wcoefs[2] - vcoefs[3])
    curl.append(ucoefs[3] - wcoefs[1])
    curl.append(vcoefs[1] - ucoefs[2])
    return curl


# data = getRectangularGridWithVectors(10, 10, 10)
# test_cell = data[5, 5, 5]
# solve(test_cell)
# print(str(test_cell.polyfitAverage[0]))


