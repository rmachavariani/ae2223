import gridConstructionSphereFast as gr
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import EnsemblePolyfit as ens
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# parameters
nrOfParticles = 50000  # None

pitch = [20]  # [10,15,20]
radius = [20]  # [10,15,20]
types = ["nor"]  # ,"gauss","poly","dog"
planes = ["yz", "xz", "xy"]

# load the grids
grids = gr.allgrid(pitch, radius, nrOfParticles)

t1 = time.time()

# loop over grids to calculate averages
for grid in grids:

    # loop over current grid
    for i in range(np.size(grid, axis=0)):
        for j in range(np.size(grid, axis=1)):
            for k in range(np.size(grid, axis=2)):

                # current bin
                thisBin = grid[i, j, k]

                # calculate normal average
                if types.count("nor") == 1:
                    thisBin.calculateNormalAverage()

                # calculate variance and deviation
                if types.count("gauss") == 1:
                    thisBin.calculateStandardDeviation()
                    thisBin.calculateVariance()

                    # calculate gaussian average
                    thisBin.calculateGaussianAverage()

                # calculate dog gaussian average

                # calculate poly fit
                if len(thisBin.vectors) > 20 and types.count("poly") == 1:
                    ens.solve(thisBin)

t2 = time.time()

print()
print("Total time for taking averages: ", round(t2 - t1, 2))
print()
print("----------------------------------------")


def plotting(planes, typ, location, unit, number):
    t1 = time.time()
    for grid in grids:
        # figure, axes = plt.subplots(nrows=1, ncols=3)

        levels = np.linspace(-10, 25, number)

        uyz_values, vyz_values, wyz_values = [], [], []
        xyz_values, yyz_values, zyz_values = [], [], []
        uxz_values, vxz_values, wxz_values = [], [], []
        xxz_values, yxz_values, zxz_values = [], [], []
        uxy_values, usxy, vxy_values, vsxy, wxy_values = [], [], [], [], []
        xxy_values, yxy_values, zxy_values = [], [], []

        if planes.count("yz"):
            for i in range(np.size(grid, axis=0)):
                if grid[i, 0, 0].x >= location:
                    a = i
                    break
            for i in range(np.size(grid, axis=1)):
                for j in range(np.size(grid, axis=2)):
                    thisBin = grid[a, i, j]
                    xyz_values.append(thisBin.x * 0.001)
                    yyz_values.append(thisBin.y * 0.001)
                    zyz_values.append(thisBin.z * 0.001)

                    if typ == "poly":
                        lst = thisBin.polyfitAverage
                    elif typ == "nor":
                        lst = thisBin.averageNormal
                    elif typ == "gauss":
                        lst = thisBin.averageGauss
                    elif typ == "dog":
                        lst = thisBin.dogGaussianAverage
                    elif typ == "vor":
                        lst = thisBin.vorticity
                    if lst != []:
                        uyz_values.append(lst[0])
                        vyz_values.append(lst[1])
                        wyz_values.append(lst[2])
                    else:
                        uyz_values.append(None)
                        vyz_values.append(None)
                        wyz_values.append(None)
        if planes.count("xz") == 1:
            for i in range(np.size(grid, axis=1)):
                if grid[0, i, 0].y >= location:
                    b = i
                    break
            for i in range(np.size(grid, axis=0)):
                for j in range(np.size(grid, axis=2)):
                    thisBin = grid[i, b, j]
                    xxz_values.append(thisBin.x * 0.001)
                    yxz_values.append(thisBin.y * 0.001)
                    zxz_values.append(thisBin.z * 0.001)

                    if typ == "poly":
                        lst = thisBin.polyfitAverage
                    elif typ == "nor":
                        lst = thisBin.averageNormal
                    elif typ == "gauss":
                        lst = thisBin.averageGauss
                    elif typ == "dog":
                        lst = thisBin.dogGaussianAverage
                    elif typ == "vor":
                        lst = thisBin.vorticity

                    if lst != []:
                        uxz_values.append(lst[0])
                        vxz_values.append(lst[1])
                        wxz_values.append(lst[2])
                    else:
                        uxz_values.append(None)
                        vxz_values.append(None)
                        wxz_values.append(None)
        if planes.count("xy") == 1:
            for i in range(np.size(grid, axis=2)):
                if grid[0, 0, i].z >= location:
                    c = i
                    break
            for i in range(np.size(grid, axis=0)):
                for j in range(np.size(grid, axis=1)):
                    thisBin = grid[i, j, c]
                    xxy_values.append(thisBin.x * 0.001)
                    yxy_values.append(thisBin.y * 0.001)
                    zxy_values.append(thisBin.z * 0.001)

                    if typ == "poly":
                        lst = thisBin.polyfitAverage
                    elif typ == "nor":
                        lst = thisBin.averageNormal
                    elif typ == "gauss":
                        lst = thisBin.averageGauss
                    elif typ == "dog":
                        lst = thisBin.dogGaussianAverage
                    elif typ == "vor":
                        lst = thisBin.vorticity

                    if lst != []:
                        uxy_values.append(lst[0])
                        usxy.append(lst[0])
                        vxy_values.append(lst[1])
                        vsxy.append(lst[1])
                        wxy_values.append(lst[2])
                    else:
                        uxy_values.append(None)
                        usxy.append(0)
                        vxy_values.append(None)
                        vsxy.append(0)
                        wxy_values.append(None)

        figure, axes = plt.subplots(nrows=2, ncols=2, gridspec_kw={
            'width_ratios': [1, ((max(yyz_values) - min(yyz_values)) / (max(xxy_values) - min(xxy_values)))],
            'height_ratios': [1, ((max(yxy_values) - min(yxy_values)) / (max(zxz_values) - min(zxz_values)))]})
        if planes.count("yz") == 1:
            ayz_lst = np.array(xyz_values).reshape(np.size(grid, axis=1), np.size(grid, axis=2))
            byz_lst = np.array(yyz_values).reshape(np.size(grid, axis=1), np.size(grid, axis=2))
            cyz_lst = np.array(zyz_values).reshape(np.size(grid, axis=1), np.size(grid, axis=2))
            dyz_lst = np.array(uyz_values).reshape(np.size(grid, axis=1), np.size(grid, axis=2))
            loc = " " + str(round(grid[a, 0, 0].x * 0.001, 4))
            plot = axes[0][1].contourf(byz_lst, cyz_lst, dyz_lst, levels=levels, cmap=cm.RdBu_r)
            axes[0][1].set_title("yz plane at x=" + loc + "[m]")
            axes[0][1].set_xlabel("y [m]")
            axes[0][1].set_ylabel("z [m]")
            axes[0][1].grid()
        if planes.count("xz") == 1:
            axz_lst = np.array(xxz_values).reshape(np.size(grid, axis=0), np.size(grid, axis=2))
            bxz_lst = np.array(yxz_values).reshape(np.size(grid, axis=0), np.size(grid, axis=2))
            cxz_lst = np.array(zxz_values).reshape(np.size(grid, axis=0), np.size(grid, axis=2))
            dxz_lst = np.array(uxz_values).reshape(np.size(grid, axis=0), np.size(grid, axis=2))
            loc = " " + str(round(grid[0, b, 0].y * 0.001, 4))
            plot = axes[0][0].contourf(axz_lst, cxz_lst, dxz_lst, levels=levels, cmap=cm.RdBu_r)
            axes[0][0].set_title("xz plane at y=" + loc + "[m]")
            axes[0][0].set_xlabel("x [m]")
            axes[0][0].set_ylabel("z [m]")
            axes[0][0].grid()
        if planes.count("xy") == 1:
            axy_lst = np.array(xxy_values).reshape(np.size(grid, axis=0), np.size(grid, axis=1))
            bxy_lst = np.array(yxy_values).reshape(np.size(grid, axis=0), np.size(grid, axis=1))
            cxy_lst = np.array(zxy_values).reshape(np.size(grid, axis=0), np.size(grid, axis=1))
            dxy_lst = np.array(uxy_values).reshape(np.size(grid, axis=0), np.size(grid, axis=1))
            loc = " " + str(round(grid[0, 0, c].z * 0.001, 4))
            plot = axes[1][0].contourf(axy_lst, bxy_lst, dxy_lst, levels=levels, cmap=cm.RdBu_r)
            axes[1][0].set_title("xy plane at z=" + loc + "[m]")
            axes[1][0].set_xlabel("x [m]")
            axes[1][0].set_ylabel("y [m]")
            axes[1][0].quiver(xxy_values, yxy_values, usxy, vsxy, scale=280, color='Black')
            axes[1][0].grid()
        axes[1][1].axis('off')

        # plt.autoscale(enable=True, axis='both', tight=None)
        t2 = time.time()
        print("Time to order values", round(t2 - t1, 3), "sec")

        names = ["Top hat filter ensemble averaging", "Gaussian ensemble averaging",
                 "Polynomial fit ensemble averaging"]
        if typ == "nor":
            name = names[0]
        elif typ == "gauss":
            name = names[1]
        elif typ == "poly":
            name = names[2]
        cb = figure.colorbar(plot)
        cb.set_label(unit)
        figure.suptitle("Velocity heatmap using " + name)
        # axes.title(typ+" Velocity Heatmap in "+plane+loc+"[mm]")
        plt.show()

        # fig = plt.figure(figsize=(10,6))
        # ax1 = fig.add_subplot(111, projection='3d')

        # ax1.set_title('2D contour plots in 3D plot for'+plane)
        # surf1 = ax1.plot_surface(a_lst, b_lst, c_lst, facecolors=cm.jet(d_lst))
        # fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=5)
        # ax1.set_xlabel('x [mm]')
        # ax1.set_ylabel('y [mm]')
        # ax1.set_zlabel('z [mm]')
        # plt.show()


for typ in types:
    plotting(planes, typ, 5, "u [m/s]", 10)

# plotting("yz","poly",0)
# plotting("xz","poly",0)
# plotting("xy","poly",10)
# plotting("yz","vor",0)
