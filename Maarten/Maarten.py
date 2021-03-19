from gridConstruction import *
import numpy as np


class Polyfit:

    # A polyomial fit on a velocity component in a 3 dimensional domain
    def __init__(self,data):
        self.vectors = data

    def basis(self,x,y,z):
        v_comp = np.array([1,x,y,z,x*y,x*z,y*z,x**2,y**2,z**2])
        return v_comp

    def design_matrix(self):
        i = 0
        matrix = np.empty(np.size(self.vectors), dtype=object)
        for vector in self.vectors:
            matrix[i] = self.basis(vector[0],vector[1],vector[2])
            i += 1
        return matrix

    def least_squares(self,design,data):




test_bin = getRectangularGridWithVectors(10, 10, 10)[5, 5, 5].vectors
particles_val = np.empty((len(test_bin), 6))
for i in range(len(test_bin)):
    x = test_bin[i].x
    y = test_bin[i].y
    z = test_bin[i].z
    u = test_bin[i].u
    v = test_bin[i].v
    w = test_bin[i].w

    particles_val[i] = [x, y, z, u, v, w]


testObject = Polyfit(particles_val[:,0:4])
test = testObject.design_matrix()
print(test)