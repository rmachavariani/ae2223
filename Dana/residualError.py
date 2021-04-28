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
    bin = grid[41,23,28]

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
    y_local = []
    z_local = []
    for x_value, y_value, z_value in zip(x_values,y_values,z_values):
        x_loc = (x_value - bin.x)
        y_loc = (y_value - bin.y)
        z_loc = (z_value - bin.z)

        x_local.append(x_loc)
        y_local.append(y_loc)
        z_local.append(z_loc)


    #normal average plot
    bin.calculateNormalAverage()
    normal_average_u = bin.normalAverage[0][0]
    normal_average_u = np.array([normal_average_u for i in range(len(x_values))])
    plt.plot(x_local, normal_average_u, label='Normal Average')

    #gaussian average plot
    bin.calculateStandardDeviation()
    bin.calculateVariance()
    bin.calculateGaussianAverage()
    gaussian_average_u = bin.gaussianAverage[0][0]
    gaussian_average_u = np.array([gaussian_average_u for i in range(len(x_values))])
    plt.plot(x_local, gaussian_average_u, label='Gaussian Average')

    #polynomial fit plot
    if len(bin.vectors) > 10:
        ens.solve(bin)
        poly_u = bin.fitU
        plt.plot(x_local, [poly_u(dx=x,dy=0,dz=0) for x,y,z in zip(x_local, y_local, z_local)], 'b+', label='Polyfit')

    # x_local = [coord * 1000 for coord in x_local]
    plt.scatter(x_local, u_values, s=10)
    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.axis('square')
    plt.ylim(top=20)
    plt.legend()
    plt.show()


