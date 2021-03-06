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