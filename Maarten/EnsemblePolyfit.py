from gridConstruction import *
from PolyClass import *
import numpy as np






test_bin = getRectangularGridWithVectors(10, 10, 10)[5, 5, 5]

testpolyfit = poly(test_bin)
print(testpolyfit.coeffs())





# def polyfit_ensemble(gridsize):
#     dataset = getRectangularGridWithVectors(gridsize[0],gridsize[1],gridsize[2])
#     for bin in dataset:
#         vector_array =