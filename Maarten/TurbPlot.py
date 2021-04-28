import gridConstructionSphereFast as gr
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import EnsemblePolyfit as ens
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
matplotlib.use("Qt5Agg")

# parameters
nrOfParticles = None

pitch = [30]#[10,15,20]
radius = [30]#[10,15,20]

# load the grids
grids = gr.allgrid(pitch,radius,nrOfParticles)


t1 = time.time()

# loop over grids to calculate averages
for grid in grids:

    # loop over current grid
    for i in range(np.size(grid,axis=0)):
        for j in range(np.size(grid,axis=1)):
            for k in range(np.size(grid,axis=2)):

                # current bin
                thisBin = grid[i,j,k]

                # calculate normal average
                #thisBin.calculateNormalAverage()

                # calculate variance and deviation
                #thisBin.calculateStandardDeviation()
                #thisBin.calculateVariance()

                # calculate gaussian average
                #thisBin.calculateGaussianAverage()

                # calculate dog gaussian average


                # calculate poly fit
                if len(thisBin.vectors) > 50:
                    ens.solve(thisBin)

t2 = time.time()

print()
print("Total time for taking averages: ", round(t2-t1,2))
print()
print("----------------------------------------")


def plotting(plane,location, locstr, number):
    t1 = time.time()
    for grid in grids:
        t_values= []
        x_values, y_values, z_values = [],[],[]

        levels = np.linspace(0, 40, number)

        if plane == "yz":
            for i in range(np.size(grid,axis=0)):
                if grid[i,0,0].x >=location:
                    a = i
                    break
            print(grid[a,0,0].x)
            for i in range(np.size(grid,axis=1)):
                for j in range(np.size(grid,axis=2)):
                    thisBin = grid[a,i,j]
                    x_values.append(thisBin.x)
                    y_values.append(thisBin.y)
                    z_values.append(thisBin.z)

                    if thisBin.polyfitAverage != []:
                        t_values.append(thisBin.turb_eng)
                    else:
                        t_values.append(0)
        elif plane == "xz":
            for i in range(np.size(grid,axis=1)):
                if grid[0,i,0].y >=location:
                    a = i
                    break
            print(grid[0,a,0].y)
            for i in range(np.size(grid,axis=0)):
                for j in range(np.size(grid,axis=2)):
                    thisBin = grid[i,a,j]
                    x_values.append(thisBin.x)
                    y_values.append(thisBin.y)
                    z_values.append(thisBin.z)

                    if thisBin.polyfitAverage != []:
                        t_values.append(thisBin.turb_eng)
                    else:
                        t_values.append(0)
        elif plane == "xy":
            for i in range(np.size(grid,axis=2)):
                if grid[0,0,i].z >=location:
                    a = i
                    break
            print(grid[0,0,a].z)
            for i in range(np.size(grid,axis=0)):
                for j in range(np.size(grid,axis=1)):
                    thisBin = grid[i,j,a]
                    x_values.append(thisBin.x)
                    y_values.append(thisBin.y)
                    z_values.append(thisBin.z)

                    if thisBin.polyfitAverage != []:
                        t_values.append(thisBin.turb_eng)
                    else:
                        t_values.append(0)

        if plane == "yz":
            b_lst = np.array(y_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))
            c_lst = np.array(z_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))
            e_lst = np.array(t_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))
            plt.contourf(b_lst, c_lst, e_lst, 100, levels=levels, cmap=cm.RdBu_r)
        elif plane == "xz":
            a_lst = np.array(x_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))
            c_lst = np.array(z_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))
            e_lst = np.array(t_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))
            plt.contourf(a_lst, c_lst, e_lst, 100, levels=levels, cmap=cm.RdBu_r)
        elif plane == "xy":
            a_lst = np.array(x_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))
            b_lst = np.array(y_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))
            e_lst = np.array(t_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))
            plt.contourf(a_lst, b_lst, e_lst, 100, levels=levels, cmap=cm.RdBu_r)

        t2 = time.time()
        print("Time to order values",round(t2-t1,3),"sec")

        cb = plt.colorbar()
        plt.axis("equal")
        plt.grid()
        cb.set_label("K [J/kg]")
        plt.title("Turbulent Kinetic Energy "+ plane + " at " + locstr + str(location))
        plt.show()


plotting("yz",0, " x = ", 9)
plotting("xz",0, " y = ", 9)
plotting("xy",0, " z = ", 9)

