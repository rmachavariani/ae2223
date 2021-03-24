import numpy as np
class poly:

    # A polyomial fit on a velocity component in a 3 dimensional domain
    def __init__(self,bin):
        self.bin = bin
        self.vectors = bin.vectors
        self.designmatrix = np.empty((len(self.vectors),10))
        self.particles_data = np.empty((len(self.vectors), 6))
        self.ucoeffs = None
        self.vcoeffs = None
        self.wcoeffs = None

    def fetch_vectors(self):
        if self.vectors == [] or () or None:
            print("No particles")
        else:
            for i in range(len(self.vectors)):
                x = self.vectors[i].x
                y = self.vectors[i].y
                z = self.vectors[i].z
                u = self.vectors[i].u
                v = self.vectors[i].v
                w = self.vectors[i].w

                self.particles_data[i] = [x, y, z, u, v, w]

    def basis(self,x,y,z):
        x = (x - self.bin.x) / 1000
        y = (y - self.bin.y) / 1000
        z = (x - self.bin.z) / 1000
        v_comp = [1.,x ,y ,z ,x*y ,x*z ,y*z ,x**2 ,y**2 ,z**2]
        return v_comp

    def design_matrix(self):
        for i in range(np.size(self.particles_data, axis=0)):
            self.designmatrix[i] = self.basis(self.particles_data[0],self.particles_data[1],self.particles_data[2])


    def coeffs(self):
        self.fetch_vectors()
        self.design_matrix()
        c = 3
        for i in [self.ucoeffs, self.vcoeffs, self.wcoeffs]:
            i = np.linalg.lstsq(self.designmatrix, self.vectors[:,c])
            c += 1
        return self.ucoeffs, self.vcoeffs, self.wcoeffs