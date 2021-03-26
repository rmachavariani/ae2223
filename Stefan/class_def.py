import numpy as np
import matplotlib.pyplot as plt

class vector:

    def __init__(self, particleData):
        self.x = particleData[0]
        self.y = particleData[1]
        self.z = particleData[2]
        self.u = particleData[3]
        self.v = particleData[4]
        self.w = particleData[5]


class gridBin:
    # bin radius and number of bins in each direction
    # for spherical bin (static member belonging to class)
    radius = 0

    # number of bins
    nrBinsX = 0
    nrBinsY = 0
    nrBinsZ = 0

    # bins width for rectangular bins (static class members)
    widthX = 0
    widthY = 0
    widthZ = 0

    xMin = 0
    xMax = 0
    yMin = 0
    yMax = 0
    zMin = 0
    zMax = 0

    def __init__(self, x, y, z):

        # coordinate of center of bin
        self.x = x
        self.y = y
        self.z = z

        # list of vectors belonging to the bin
        self.vectors = []
        self.NormalAverage = []
        self.GaussianAverage = []

    def addVector(self, vector):

        self.vectors.append(vector)

    def calculateNormalAverage(self):
        self.averageU = 0
        self.averageV = 0
        self.averageW = 0
        sumU = 0
        sumV = 0
        sumW = 0

        if np.size(self.vectors) >= 1:
            for vector in self.vectors:
                sumU += vector.u
                sumV += vector.v
                sumW += vector.w

            nrVectors = np.size(self.vectors)

            self.averageU = sumU / nrVectors
            self.averageV = sumV / nrVectors
            self.averageW = sumW / nrVectors

        self.NormalAverage.append([self.averageU, self.averageV, self.averageW])

    def calculateStandardDeviation(self):

        diffU = 0
        diffV = 0
        diffW = 0

        # there is only one vector so that would mean standard deviation = 0
        if np.size(self.vectors) >= 1:
            for vector in self.vectors:
                diffU += ((vector.u - self.averageU) ** 2)
                diffV += ((vector.v - self.averageV) ** 2)
                diffW += ((vector.w - self.averageW) ** 2)

            nrVectors = np.size(self.vectors)

            self.sdU = np.sqrt(diffU / nrVectors)
            self.sdV = np.sqrt(diffV / nrVectors)
            self.sdW = np.sqrt(diffW / nrVectors)

    def calculateVariance(self):
        self.varU = 0
        self.varV = 0
        self.varW = 0

        if np.size(self.vectors) >= 1:
            #            for vector in self.vectors:
            self.varU += self.sdU ** 2
            self.varV += self.sdV ** 2
            self.varW += self.sdW ** 2

        return self.varU, self.varV, self.varW

    def calculateGaussianAverage(self):
        gaussSumU, gaussSumV, gaussSumW = 0, 0, 0
        gaussianWeightedU, gaussianWeightedV, gaussianWeightedW = 0, 0, 0

        if np.size(self.vectors) > 1:
            # use numpy to calculate the mean and standard deviation
            for vector in self.vectors:
                # calculate for every velocity component the normal value and sum this together for the bin
                self.weightU = (1 / (self.sdU * np.sqrt(2 * np.pi))) * np.exp(
                    (-1 / 2) * (((vector.u - self.averageU) / self.sdU) ** 2))
                # multiply the gaussian weight with the velocity component and sum all in the bin
                gaussianWeightedU += self.weightU * vector.u
                # sum all gaussian weights within the bin
                gaussSumU += self.weightU

                self.weightV = (1 / (self.sdV * np.sqrt(2 * np.pi))) * np.exp(
                    (-1 / 2) * (((vector.v - self.averageV) / self.sdV) ** 2))
                gaussianWeightedV += self.weightV * vector.v
                gaussSumV += self.weightV

                self.weightW = (1 / (self.sdW * np.sqrt(2 * np.pi))) * np.exp(
                    (-1 / 2) * (((vector.w - self.averageW) / self.sdW) ** 2))
                gaussianWeightedW += self.weightW * vector.w
                gaussSumW += self.weightW

            # devide the sum of the normal*velocity by the sum of the normal values to get the average gaussian velocity
            self.gaussU = gaussianWeightedU / (gaussSumU)
            self.gaussV = gaussianWeightedV / (gaussSumV)
            self.gaussW = gaussianWeightedW / (gaussSumW)

        elif np.size(self.vectors) == 1:
            # there is only one vector so that would mean standard deviation = 0 and Gaussian
            # method can not be appolied so the average velocity is the velocity component of the particle
            for vector in self.vectors:
                self.gaussU = vector.u
                self.gaussV = vector.v
                self.gaussW = vector.w

        elif np.size(self.vectors) == 0:
            # for an empty bin the velocity component is just 0
            self.gaussU = 0
            self.gaussV = 0
            self.gaussW = 0

        self.GaussianAverage.append([self.gaussU, self.gaussV, self.gaussW])


    def calculateNewGaussianAverage(self,datavarU, datavarV, datavarW, bigdatavarU, bigdatavarV, bigdatavarW):
        gaussSumU, gaussSumV, gaussSumW = 0, 0, 0
        gaussianWeightedU, gaussianWeightedV, gaussianWeightedW = 0, 0, 0

        if np.size(self.vectors) > 1:
            # use numpy to calculate the mean and standard deviation
            for vector in self.vectors:
                # calculate for every velocity component the normal value and sum this together for the bin
                weightU = self.weightU
                weightU2 = (1 / (np.sqrt(bigdatavarU) * np.sqrt(2 * np.pi))) * np.exp(-((vector.u - self.averageU) ** 2) / (2 * bigdatavarU))
                weightUtot = weightU2 - weightU

                # multiply the gaussian weight with the velocity component and sum all in the bin
                gaussianWeightedU += weightUtot * vector.u
                # sum all gaussian weights within the bin
                gaussSumU += weightUtot

                # calculate for every velocity component the normal value and sum this together for the bin
                weightV = self.weightV
                weightV2 = (1 / (np.sqrt(bigdatavarV) * np.sqrt(2 * np.pi))) * np.exp(-((vector.v - self.averageV) ** 2) / (2 * bigdatavarV))
                weightVtot = weightV2 - weightV

                # multiply the gaussian weight with the velocity component and sum all in the bin
                gaussianWeightedV += weightVtot * vector.v
                # sum all gaussian weights within the bin
                gaussSumV += weightVtot

                # calculate for every velocity component the normal value and sum this together for the bin
                weightW = self.weightW
                weightW2 = (1 / (np.sqrt(bigdatavarW) * np.sqrt(2 * np.pi))) * np.exp(-((vector.w - self.averageW) ** 2) / (2 * bigdatavarW))
                weightWtot = weightW2 - weightW

                # multiply the gaussian weight with the velocity component and sum all in the bin
                gaussianWeightedW += weightWtot * vector.w
                # sum all gaussian weights within the bin
                gaussSumW += weightWtot

            # devide the sum of the normal*velocity by the sum of the normal values to get the average gaussian velocity
            self.newgaussU = gaussianWeightedU / (gaussSumU)
            self.newgaussV = gaussianWeightedV / (gaussSumV)
            self.newgaussW = gaussianWeightedW / (gaussSumW)

        elif np.size(self.vectors) == 1:
            # there is only one vector so that would mean standard deviation = 0 and Gaussian
            # method can not be appolied so the average velocity is the velocity component of the particle
            for vector in self.vectors:
                self.newgaussU = vector.u
                self.newgaussV = vector.v
                self.newgaussW = vector.w

        elif np.size(self.vectors) == 0:
            # for an empty bin the velocity component is just 0
            self.newgaussU = 0
            self.newgaussV = 0
            self.newgaussW = 0

        self.newGaussianAverage.append([self.newgaussU, self.newgaussV, self.newgaussW])
