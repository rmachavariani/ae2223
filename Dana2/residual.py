import matplotlib

# import gridConstructionSphereFast as gr
import gridConstructionSphereFast_Big as gr

import numpy as np
import EnsemblePolyfit as ens
import matplotlib.pyplot as plt

del matplotlib.font_manager.weight_dict['roman']
# matplotlib.font_manager._rebuild()


# parameters
nrOfParticles = None

pitch = [15]
radius = [7.5]
# load the grids
grids = gr.allgrid(pitch,radius,nrOfParticles)


# def get3DPlot(x_positions,y_positions,z_positions):
#     fig = plt.figure()
#     ax = mplot3d.Axes3D(fig)
#
#     # Load the STL files and add the vectors to the plot
#     car_mirror = mesh.Mesh.from_file('Dana2/modelScaled.stl')
#     car_mirror.rotate([0, 0, 1], math.radians(180))
#     ax.add_collection3d(mplot3d.art3d.Poly3DCollection(car_mirror.vectors,edgecolor='k',facecolors='w', linewidths=1, alpha=0.5))
#     # Auto scale to the mesh size
#     scale = car_mirror.points.flatten()
#     ax.auto_scale_xyz(scale, scale, scale)
#
#     x_scaled = []
#     y_scaled = []
#     z_scaled = []
#
#     for x,y,z in zip(x_positions, y_positions, z_positions):
#         x_scaled.append(x/1000)
#         y_scaled.append(y/1000)
#         z_scaled.append(z/1000)
#
#     ax.set_xlabel('X [m]')
#     ax.set_ylabel('Y [m]')
#     ax.set_zlabel('Z [m]')
#     ax.scatter(x_scaled[0::100], y_scaled[0::100], z_scaled[0::100], s=5)
#     plt.show()


def getAveragingPlot(bin,i,j,k):
    vectors = bin.vectors
    fig = plt.figure()

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
    for x_value, y_value, z_value in zip(x_values, y_values, z_values):
        x_loc = (x_value - bin.x)
        y_loc = (y_value - bin.y)
        z_loc = (z_value - bin.z)

        x_local.append(x_loc)
        y_local.append(y_loc)
        z_local.append(z_loc)

    y_local_over_H = []
    for y_local_mm in y_local:
        y_local_over_H.append(y_local_mm/0.75)


    # normal average plot
    bin.calculateNormalAverage()
    normal_average_u = bin.normalAverage[0][0]
    normal_average_u = np.array([normal_average_u for i in range(len(y_local_over_H))])
    plt.plot(y_local_over_H, normal_average_u, label='Top Hat Average', color='y')

    # gaussian average plot
    bin.calculateStandardDeviation()
    bin.calculateVariance()
    bin.calculateGaussianAverage()
    gaussian_average_u = bin.gaussianAverage[0][0]
    gaussian_average_u = np.array([gaussian_average_u for i in range(len(y_local_over_H))])
    plt.plot(y_local_over_H, gaussian_average_u, label='Gaussian Average', color='m')

    # polynomial fit plot
    if len(bin.vectors) > 20:
        ens.solve(bin)
        poly_u = bin.fitU
        plt.plot(y_local_over_H, [poly_u(dx=0, dy=y, dz=0) for x, y, z in zip(x_local, y_local, z_local)], 'b+',
                 label='Polyfit')

    plt.scatter(y_local_over_H, u_values, s=5, color='c')

    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.7)

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

    plt.axvline(linewidth=1, color='r', linestyle='dashed')

    # plt.ylim(top=20)
    plt.xlim(-1,1)
    # plt.axis('equal')
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = 16
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.xlabel('y/H')
    plt.ylabel('Velocity [m/s]')
    file_name = 'plots3/' + str(i) + '_' + str(j) + '_' + str(k) + '.png'
    title = 'Bin ' + str(bin.x) + ',' + str(bin.y) + ',' + str(bin.z)
    plt.title(title)
    # plt.legend(loc="upper left",bbox_to_anchor=(1.04,1))
    # plt.tight_layout()
    plt.savefig(file_name)
    fig.clear()
    plt.close(fig)


for grid in grids:
    x_positions = []
    y_positions = []
    z_positions = []

    for i in range(0,np.size(grid,axis=0)):
        for j in range(0,np.size(grid,axis=1)):
            for k in range(0,np.size(grid,axis=2)):

                bin = grid[i,j,k]

                vectors = bin.vectors

                getAveragingPlot(bin,i,j,k)

                # #current bin
                # bins = [grid[39,23,23], grid[59,38,14], grid[70,38,14]]
                # getAveragingPlot(bins)



    # x_positions.append(bin.x)
    # y_positions.append(bin.y)
    # z_positions.append(bin.z)


    # get3DPlot(x_positions,y_positions,z_positions)

