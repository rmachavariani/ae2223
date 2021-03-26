import functools.partial as part
import numpy as np

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
