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

    def addVector(self, vector):

        self.vectors.append(vector)

    def calculateNormalAverage(self):

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

    def calculateGaussianAverage(self):
        gaussSumU, gaussSumV, gaussSumW = 0, 0, 0
        gaussianWeightedU, gaussianWeightedV, gaussianWeightedW = 0, 0, 0

        if np.size(self.vectors) > 1:
            # use numpy to calculate the mean and standard deviation
            for vector in self.vectors:
                # calculate for every velocity component the normal value and sum this together for the bin
                weightU = (1 / (self.sdU * np.sqrt(2 * np.pi))) * np.exp(
                    (-1 / 2) * (((vector.u - self.averageU) / self.sdU) ** 2))
                # multiply the gaussian weight with the velocity component and sum all in the bin
                gaussianWeightedU += weightU * vector.u
                # sum all gaussian weights within the bin
                gaussSumU += weightU

                weightV = (1 / (self.sdV * np.sqrt(2 * np.pi))) * np.exp(
                    (-1 / 2) * (((vector.v - self.averageV) / self.sdV) ** 2))
                gaussianWeightedV += weightV * vector.v
                gaussSumV += weightV

                weightW = (1 / (self.sdW * np.sqrt(2 * np.pi))) * np.exp(
                    (-1 / 2) * (((vector.w - self.averageW) / self.sdW) ** 2))
                gaussianWeightedW += weightW * vector.w
                gaussSumW += weightW

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
    def histogram(self):
        lst_u = []
        lst_v = []
        lst_w = []
        for vector in self.vectors:
            lst_u.append(vector.u)
            lst_v.append(vector.v)
            lst_w.append(vector.w)

        plt.figure()
        #d_v = np.linspace(min(min(lst_u),min(lst_v),min(lst_w))+10, max(max(lst_u),max(lst_v),max(lst_w))-10, 30, dtype=np.float64)
        plt.subplot(1,1,1)
        plt.title("Histogram bin:"+str(5)+","+str(5)+","+str(5))
        plt.hist([lst_u,lst_v,lst_w], bins = 30, label=["u","v","w"])
        plt.legend(loc='upper right')
        plt.xlabel("Velocity [m/s]")
        plt.ylabel("Count")
        plt.show()