import Dana.gridConstructionSphereFast as gr
import numpy as np
import time
import matplotlib.pyplot as plt
import Dana.EnsemblePolyfit as ens

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
    y_values = []
    z_values = []


    for vector in vectors:
        u_values.append(vector.u)
        x_values.append(vector.x)
        y_values.append(vector.y)
        z_values.append(vector.z)

    x_local = []
    for value in range(len(x_values)):
        x_loc = (value - bin.x)
        x_local.append(x_loc)


    # #normal average plot
    # bin.calculateNormalAverage()
    # normal_average_u = bin.normalAverage[0][0]
    # normal_average_u = np.array([normal_average_u for i in range(len(x_values))])
    # plt.plot(x_values, normal_average_u)
    #
    # #gaussian average plot
    # bin.calculateStandardDeviation()
    # bin.calculateVariance()
    # bin.calculateGaussianAverage()
    # gaussian_average_u = bin.gaussianAverage[0][0]
    # gaussian_average_u = np.array([gaussian_average_u for i in range(len(x_values))])
    # plt.plot(x_values, gaussian_average_u)

    #polynomial fit plot
    if len(bin.vectors) > 10:
        ens.solve(bin)
        poly_u = bin.fitU
        print(poly_u)
        plt.plot(x_local, [poly_u(dx=item,dy=0,dz=0) for item in x_local], 'b+', label='Polyfit')

    plt.scatter(x_local, u_values)
    plt.show()


