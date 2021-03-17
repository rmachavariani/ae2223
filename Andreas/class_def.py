
class vector:
    
    def __init__(self,particleData):
        self.x = particleData[0]
        self.y = particleData[1]
        self.z = particleData[2]
        self.u = particleData[3]
        self.v = particleData[4]
        self.w = particleData[5]
        
class gridBin:
    

    # bins width (static class members)
    widthX = 0
    widthY = 0
    widthZ = 0
    
    def __init__(self,x,y,z):

        # coordinate of center of bin
        self.x = x
        self.y = y
        self.z = z

        # list of vectors belonging to the bin
        self.vectors = []
        
    def addVector(self,vector):
        
        self.vectors.append(vector)