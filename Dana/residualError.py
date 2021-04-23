import Dana.gridConstructionSphereFast as gr
import numpy as np
import time
import matplotlib.pyplot as plt

# parameters
nrOfParticles = None

pitch = [10]
radius = [10]
# load the grids
grids = gr.allgrid(pitch,radius,nrOfParticles)

t1 = time.time()
for grid in grids:

    # current bin
    bin = grid[41,20,28]

    vectors = bin.vectors
    u_values = []
    x_values = []


    for vector in vectors:
        u_values.append(vector.u)
        x_values.append(vector.x)

    mymodel = np.poly1d(np.polyfit(x_values, u_values, 3))

    myline = np.linspace(1, 22, 100)

    plt.plot(myline, mymodel(myline))

    plt.scatter(x_values, x_values)
    plt.show()


