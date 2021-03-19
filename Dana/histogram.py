import matplotlib.pyplot as plt
import numpy as np
from Dana.class_def import *
from Dana.gridConstruction import getRectangularGridWithVectors

def visualize():
    data = getRectangularGridWithVectors(10,10,10)

    sizeX = np.size(data,axis=0)
    sizeY = np.size(data,axis=1)
    sizeZ = np.size(data,axis=2)

    lenArray = []


    print(lenArray)
    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                lenArray.append([i, j, k, len(data[i][j][k].vectors)])


    lenArray = np.array(lenArray)
    largestKeys = np.argpartition(lenArray[:,3], -4)[-4:]

    bin = lenArray[largestKeys[0]]
    bin = gridBin(bin[0], bin[1], bin[2])

    print(bin.vectors)









visualize()