from class_def import *
from gridConstruction import *
import numpy as np
#
# class Bin:
#     def __init__(self,datapoints):
#         self.vectors = np.empty(datapoints, dtype=object)
#
#
#     def addvector(self,vector,location):
#         self.vectors[location] = vector
#
#     def createdatapoints(self):
#         x = [1,1.3,1.5,1.6,1.73,1.95]
#         y = [0.47, 0.82, 0.93, 0.84, 0.75]
#         z = [1.35,1.84, 1.95,2.95,2.99]
#
#
#
#
# class Vector:
#     def __init__(self,x,y,z):
#         self.u = x**2*2.62345 + y*z*1.645783 + y**2*1.6332612345
#         self.v = y*z**2*2.62345 + y*z*1.2459 + y**2*1.6332612345
#         self.w = z*x**2*2.918374 + y*z*1.2986764 + y**2*1.6332612345

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


testcell = getRectangularGridWithVectors(10,10,10)[5,5,5][:,0:4]
testObject = Polyfit(testcell)
test = testObject.design_matrix()
print(test)