import gridConstructionSphereFast as gr
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import EnsemblePolyfit as ens
import vorticity as vort
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# parameters
nrOfParticles = None

pitch = [7.5]#[75,45,22.5,15]  # [10,6.95,26,18.081,40,27.8]#[10,20,40] 6.95,13.9,27.8
radius = [15]#[50,30,15,10]  # [5,5,13,13,20,20]#[10,15,20] 5,10,20
minParticlesForAverages = 20
types = ["poly"]  # "nor","gauss","poly","dog"
planes = ["yz", "xz", "xy"]  # "yz","xz","xy"

# load the grids
grids = gr.allgrid(pitch, radius, nrOfParticles)

t1 = time.time()

# loop over grids to calculate averages
for grid in grids:
    # loop over current grid
    for i in range(np.size(grid, axis=0)):
        nex = True
        for j in range(np.size(grid, axis=1)):
            for k in range(np.size(grid, axis=2)):
                if i % 10 ==0 and nex == True:
                    print(i)
                    nex = False
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
                if len(thisBin.vectors) > minParticlesForAverages and types.count("poly") == 1:
                    ens.solve(thisBin)

t2 = time.time()

print()
print("Total time for taking averages: ", round(t2 - t1, 2))
print()
print("----------------------------------------")


for grid in grids:

    # loop over current grid
    for i in range(np.size(grid,axis=0)):
        for j in range(np.size(grid,axis=1)):
            for k in range(np.size(grid,axis=2)):

                # current bin
                thisBin = grid[i, j, k]

                if len(thisBin.vectors) > minParticlesForAverages:
                    thisBin.vorticity = vort.vorticity(grid, thisBin, [np.size(grid, axis=0) - 1, np.size(grid, axis=1) - 1, np.size(grid, axis=2) - 1],3)


def plotting(planes, typ, location, unit, number, start, stop, title):
    t1 = time.time()
    m = 0
    n = 0
    l = 0
    for grid in grids:
        # figure, axes = plt.subplots(nrows=1, ncols=3)

        levels = np.linspace(start, stop, number)

        uyz_values, vyz_values, wyz_values = [], [], []
        xyz_values, yyz_values, zyz_values = [], [], []
        uxz_values, vxz_values, wxz_values = [], [], []
        xxz_values, yxz_values, zxz_values = [], [], []
        uxy_values, usxy, vxy_values, vsxy, wxy_values = [], [], [], [], []
        xxy_values, yxy_values, zxy_values = [], [], []

        if planes.count("yz"):
            for i in range(np.size(grid, axis=0)):
                if location == "spec":
                    if grid[i,0,0].x >=75:
                        a = i
                        break
                else:
                    if grid[i,0,0].x >=location:
                        a = i
                        break
            for i in range(np.size(grid, axis=1)):
                for j in range(np.size(grid, axis=2)):
                    thisBin = grid[a, i, j]
                    xyz_values.append(thisBin.x * 0.001)
                    yyz_values.append(thisBin.y * 0.001)
                    zyz_values.append(thisBin.z * 0.001)

                    if typ == "poly" or typ=="ke":
                        lst = thisBin.polyfitAverage
                    elif typ == "nor":
                        lst = thisBin.averageNormal
                    elif typ == "gauss":
                        lst = thisBin.averageGauss
                    elif typ == "dog":
                        lst = thisBin.dogGaussianAverage
                    elif typ == "vor":
                        lst = thisBin.vorticity

                    if typ != "ke":
                        if lst!=[] and type(lst)!=type(None):
                            if typ != "vor":
                                uyz_values.append(lst[0])
                            elif typ == "vor":
                                uyz_values.append(lst[0]*(0.15/12))
                            vyz_values.append(lst[1])
                            wyz_values.append(lst[2])
                        elif typ=="vor" and type(lst)==type(None):
                            uyz_values.append(0)
                            vyz_values.append(0)
                            wyz_values.append(0)
                        else:
                            uyz_values.append(None)
                            vyz_values.append(None)
                            wyz_values.append(None)
                    elif typ == "ke":
                        if lst != []:
                            uyz_values.append(thisBin.turb_eng)
                        else:
                            uyz_values.append(None)
        if planes.count("xz") == 1:
            for i in range(np.size(grid, axis=1)):
                if location == "spec":
                    if grid[0,i,0].y >=0:
                        b = i
                        break
                else:
                    if grid[0,i,0].y >=location:
                        b = i
                        break
            for i in range(np.size(grid, axis=0)):
                for j in range(np.size(grid, axis=2)):
                    thisBin = grid[i, b, j]
                    xxz_values.append(thisBin.x * 0.001)
                    yxz_values.append(thisBin.y * 0.001)
                    zxz_values.append(thisBin.z * 0.001)

                    if typ=="poly":
                        lst = thisBin.polyfitAverage
                    elif typ=="nor":
                        lst = thisBin.averageNormal
                    elif typ=="gauss":
                        lst = thisBin.averageGauss
                    elif typ=="dog":
                        lst = thisBin.dogGaussianAverage
                    elif typ=="vor":
                        lst = thisBin.vorticity
                    elif typ=="ke":
                        energy = thisBin.turb_eng
                        lst = thisBin.polyfitAverage

                    if typ != "ke":
                        if lst!=[] and type(lst)!=type(None):
                            if typ != "vor":
                                uxz_values.append(lst[0])
                            elif typ == "vor":
                                uxz_values.append(lst[0]*(0.15/12))
                            vxz_values.append(lst[1])
                            wxz_values.append(lst[2])
                        elif typ=="vor" and type(lst)==type(None):
                            uxz_values.append(0)
                            vxz_values.append(0)
                            wxz_values.append(0)
                        else:
                            uxz_values.append(None)
                            vxz_values.append(None)
                            wxz_values.append(None)
                    elif typ=="ke":
                        if lst != []:
                            uxz_values.append(thisBin.turb_eng)
                        else:
                            uxz_values.append(None)
        if planes.count("xy") == 1:
            for i in range(np.size(grid, axis=2)):
                if location == "spec":
                    if grid[0,0,i].z >=75:
                        c = i
                        break
                else:
                    if grid[0,0,i].z >=location:
                        c = i
                        break
            for i in range(np.size(grid, axis=0)):
                for j in range(np.size(grid, axis=1)):
                    thisBin = grid[i, j, c]
                    xxy_values.append(thisBin.x * 0.001)
                    yxy_values.append(thisBin.y * 0.001)
                    zxy_values.append(thisBin.z * 0.001)

                    if typ == "poly" or typ=="ke":
                        lst = thisBin.polyfitAverage
                    elif typ == "nor":
                        lst = thisBin.averageNormal
                    elif typ == "gauss":
                        lst = thisBin.averageGauss
                    elif typ == "dog":
                        lst = thisBin.dogGaussianAverage
                    elif typ == "vor":
                        lst = thisBin.vorticity

                    if typ != "ke":
                        if lst!=[] and type(lst)!=type(None):
                            if typ != "vor":
                                uxy_values.append(lst[0])
                            elif typ == "vor":
                                uxy_values.append(lst[0]*(0.15/12))
                            usxy.append(lst[0])
                            vxy_values.append(lst[1])
                            vsxy.append(lst[1])
                            wxy_values.append(lst[2])
                        elif typ=="vor" and type(lst)==type(None):
                            uxy_values.append(0)
                            vxy_values.append(0)
                            wxy_values.append(0)
                        elif typ=="ke" and type(lst)==type(None):
                            uxy_values.append(0)
                        else:
                            uxy_values.append(None)
                            usxy.append(0)
                            vxy_values.append(None)
                            vsxy.append(0)
                            wxy_values.append(None)
                    elif typ=="ke":
                        if lst != []:
                            uxy_values.append(thisBin.turb_eng)
                        else:
                            uxy_values.append(None)

        if typ == "nor" or typ == "gauss" or typ == "poly":
            figure, axes = plt.subplots(nrows=2, ncols=2, gridspec_kw={
                'width_ratios': [1, ((max(yyz_values) - min(yyz_values)) / (max(xxy_values) - min(xxy_values)))],
                'height_ratios': [1, ((max(yxy_values) - min(yxy_values)) / (max(zxz_values) - min(zxz_values)))]})
        elif typ == "vor" or typ=="ke":
            figure, axes = plt.subplots(nrows=1, ncols=1)

        if planes.count("yz")==1:
            ayz_lst = np.array(xyz_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))/0.15
            byz_lst = np.array(yyz_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))/0.15
            cyz_lst = np.array(zyz_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))/0.15
            if typ == "nor" or typ=="gauss" or typ == "poly":
                dyz_lst = np.array(uyz_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))
                plot = axes[0][1].contourf(byz_lst, cyz_lst, dyz_lst, levels=levels, cmap=cm.jet)
                loc = " "+str(round(grid[a,0,0].x*0.001,4))
                axes[0][1].set_xticks(np.arange(-1.0,1.1,step=0.5))
                axes[0][1].set_yticks(np.arange(0,1.6,step=0.5))
                # axes[0][1].set_xticks(np.arange(round((min(yyz_values)/0.15)*2)/2,(max(yyz_values)/0.15),step=0.5))
                # axes[0][1].set_yticks(np.arange(round((min(zyz_values)/0.15)*2)/2,(max(zyz_values)/0.15),step=0.5))
                axes[0][1].set_xlim([-1.1,1.1])
                axes[0][1].set_ylim([-0.1,1.6])
                axes[0][1].set_title("yz plane at x="+loc+"[m]",fontsize=16)
                axes[0][1].set_xlabel("y/H",fontsize=16)
                axes[0][1].set_ylabel("z/H",fontsize=16)
                axes[0][1].grid()
            elif typ == "vor" or typ == "ke":
                dyz_lst = np.array(uyz_values).reshape(np.size(grid,axis=1),np.size(grid,axis=2))
                plot = axes.contourf(byz_lst, cyz_lst, dyz_lst, levels=levels, cmap=cm.jet)
                loc = " "+str(round(grid[a,0,0].x*0.001,4))
                axes.set_xticks(np.arange(-1.0,1.1,step=0.5))
                axes.set_yticks(np.arange(0,1.6,step=0.5))
                axes.set_xlim([-1.1,1.1])
                axes.set_ylim([-0.1,1.6])
                axes.set_title("yz plane at x="+loc+"[m]",fontsize=16)
                axes.set_xlabel("y/H",fontsize=16)
                axes.set_ylabel("z/H",fontsize=16)
                axes.grid()
                
        if planes.count("xz")==1:
            axz_lst = np.array(xxz_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))/0.15
            bxz_lst = np.array(yxz_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))/0.15
            cxz_lst = np.array(zxz_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))/0.15
            if typ == "nor" or typ=="gauss" or typ == "poly":
                dxz_lst = np.array(uxz_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))
                loc = " "+str(round(grid[0,b,0].y*0.001,4))
                plot = axes[0][0].contourf(axz_lst, cxz_lst, dxz_lst, levels=levels, cmap=cm.jet)
                axes[0][0].set_xticks(np.arange(-1.5,2.1,step=0.5))
                axes[0][0].set_yticks(np.arange(0,1.6,step=0.5))
                axes[0][0].set_xlim([-1.6,2.1])
                axes[0][0].set_ylim([-0.1,1.6])
                axes[0][0].set_title("xz plane at y="+loc+"[m]",fontsize=16)
                axes[0][0].set_xlabel("x/H",fontsize=16)
                axes[0][0].set_ylabel("z/H",fontsize=16)
                axes[0][0].grid()
            elif typ == "vor" or typ == "ke":
                dxz_lst = np.array(uxz_values).reshape(np.size(grid,axis=0),np.size(grid,axis=2))
                plot = axes.contourf(axz_lst, cxz_lst, dxz_lst, levels=levels, cmap=cm.jet)
                loc = " "+str(round(grid[0,b,0].y*0.001,4))
                axes.set_xticks(np.arange(-1.5,2.1,step=0.5))
                axes.set_yticks(np.arange(0,1.6,step=0.5))
                axes.set_xlim([-1.5,2.1])
                axes.set_ylim([-0.1,1.6])
                axes.set_title("xz plane at y="+loc+"[m]",fontsize=16)
                axes.set_xlabel("x/H",fontsize=16)
                axes.set_ylabel("z/H",fontsize=16)
                axes.grid()

        if planes.count("xy")==1:
            axy_lst = np.array(xxy_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))/0.15
            bxy_lst = np.array(yxy_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))/0.15
            cxy_lst = np.array(zxy_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))/0.15
            if typ == "nor" or typ=="gauss" or typ == "poly":
                dxy_lst = np.array(uxy_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))
                loc = " "+str(round(grid[0,0,c].z*0.001,4))
                plot = axes[1][0].contourf(axy_lst, bxy_lst, dxy_lst, levels=levels, cmap=cm.jet)#cm.RdBu_r
                axes[1][0].set_title("xy plane at z="+loc+"[m]",fontsize=16)
                axes[1][0].set_xlabel("x/H",fontsize=16)
                axes[1][0].set_ylabel("y/H",fontsize=16)
                axes[1][0].quiver(axy_lst, bxy_lst, usxy, vsxy, scale=325,color='Black')
                axes[1][0].set_xticks(np.arange(-1.5,2.1,step=0.5))
                axes[1][0].set_yticks(np.arange(-1.0,1.1,step=0.5))
                axes[1][0].set_xlim([-1.6,2.1])
                axes[1][0].set_ylim([-1.1,1.1])
                axes[1][0].grid()
            elif typ == "vor" or typ == "ke":
                dxy_lst = np.array(uxy_values).reshape(np.size(grid,axis=0),np.size(grid,axis=1))
                plot = axes.contourf(axy_lst, bxy_lst, dxy_lst, levels=levels, cmap=cm.jet)
                loc = " "+str(round(grid[0,0,c].z*0.001,4))
                axes.set_xticks(np.arange(-1.5,2.1,step=0.5))
                axes.set_yticks(np.arange(-1.0,1.1,step=0.5))
                axes.set_xlim([-1.6,2.1])
                axes.set_ylim([-1.1,1.1])
                axes.set_title("xy plane at z="+loc+"[m]",fontsize=16)
                axes.set_xlabel("x/H",fontsize=16)
                axes.set_ylabel("y/H",fontsize=16)
                axes.grid()

        if typ != "vor" and typ!="ke":
            axes[1][1].axis('off')
            figure.tight_layout()

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
        elif typ == "vor":
            name = "vorticity"
        cb = figure.colorbar(plot)
        cb.set_label(unit)
        if typ == "nor" or typ == "gauss" or typ == "poly":
            figure.suptitle("Velocity heatmap using " + name)
        elif typ == "vor":
            figure.suptitle("Vorticity using " + name)

        #plt.savefig(title[l])
        plt.show()
        l+=1


filename = ['/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7/poly_50_75','/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7/poly_30_45','/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7/poly_15_225','/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7/poly_10_15']
plotting(planes, "poly", "spec", "u [m/s]", 19, -10, 25,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/trash"])

plotting("yz","ke",50,"k [J/kg]",19,0,40,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/ke_yz_50"])
plotting("xz","ke",0,"k [J/kg]",19,0,40,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/ke_xz_0"])
plotting("xy","ke",50,"k [J/kg]",19,0,40,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/ke_xy_50"])

plotting("yz","vor",30,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_yz_30"])
plotting("yz","vor",50,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_yz_50"])
plotting("yz","vor",75,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_yz_75"])
plotting("yz","vor",100,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_yz_100"])
plotting("yz","vor",150,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_yz_150"])
plotting("yz","vor",200,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_yz_200"])
plotting("xz","vor",0,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_xz_0"])
plotting("xy","vor",10,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_xy_10"])
plotting("xy","vor",75,"\u03C9(H/U)",19,-2.5,2.5,["/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Images/V7_2.7/vor_xy_75"])

# plotting("yz","vor",5,"[rad/s]",10,-1,1)
# plotting("yz","vor",15,"[rad/s]",10,-1,1)\
# plotting("yz","vor",5 ,"$\omega$*(H/U)",18,-1.5,1.5)
# plotting("yz", "vor", 10, "$\omega$*(H/U)", 18, -1.5, 1.5)
# plotting("xy", "vor", 25, "$\omega$*(H/U)", 18, -1.5, 1.5)
# plotting("xz","vor",50,"$\omega$*(H/U)",18,-1.5,1.5)


#velocity graphs


