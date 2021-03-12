import timing
#from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import math
#import earthpy as et

# Import average monthly precip to numpy array
dtype1 = np.dtype([("x",float),("y",float),("z",float),("u",float),("v",float),("w",float)])
datalocation = "/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Coding/carMirrorData.dat"
data = np.loadtxt(datalocation,dtype=dtype1)

#order the data values with x from small to large x
data_order = np.sort(data, order="x")
print("Number of particles",len(data_order))

#define outer bounds of the measurement volume
x_min, x_max = min(data["x"]), max(data["x"])
y_min, y_max = min(data["y"]), max(data["y"])
z_min, z_max = min(data["z"]), max(data["z"])

#amount of bins that are placed along the axis 
windowx, windowy, windowz = 200, 100, 75

#can make a offset were to start compared to the min and maximum location at the axis
offset = [[100,-100],[100,-100],[100,-100]]
delta_x = ((x_max+offset[0][1])-(x_min+offset[0][0]))/windowx
delta_y = ((y_max+offset[1][1])-(y_min+offset[1][0]))/windowy
delta_z = ((z_max+offset[2][1])-(z_min+offset[2][0]))/windowz

#List that define the locations of every bin (additional 1 is to include the last point as well)
x_loc = np.arange(start=x_min+offset[0][0],stop=x_max+offset[0][1]+1,step=delta_x)
y_loc = np.arange(start=y_min+offset[1][0],stop=y_max+offset[1][1]+1,step=delta_y)
z_loc = np.arange(start=z_min+offset[2][0],stop=z_max+offset[2][1]+1,step=delta_z)

#make list with all bins were the particle information can be added to
binxyz = [[[[] for k in range(windowz)] for l in range(windowy)] for i in range(windowx)] #[x:[y:[z:[particle data:()]]]]

#with the use of dividing the location of the particles by the bin width
#the rounded number results in the bin the particle has to be placed in
number_processed=0
for i in range(len(data_order)):
    #make sure that the particle is inside the measurement volume defined with offset
    if data_order[i][0]+abs(x_min)>=offset[0][0] and data_order[i][1]+abs(y_min)>=offset[1][0] and data_order[i][2]+abs(z_min)>=offset[2][0]:
        if data_order[i][0]<=x_max+offset[0][1] and data_order[i][1]<=y_max+offset[1][1] and data_order[i][2]<=z_max+offset[2][1]:
            #as there are negative coordinate points add the minimum value to get positive values
            #that can be used to determine location
            x=data_order[i][0]+abs(x_min)-offset[0][0]
            #int makes sure that the particle is assigned to the correct place in the list
            pos1 = int(x/delta_x)

            y=data_order[i][1]+abs(y_min)-offset[1][0] 
            pos2 = int(y/delta_y)

            z=data_order[i][2]+abs(z_min)-offset[2][0] 
            pos3 = int(z/delta_z)

            #the maximum coordinates will be assigned pos=window. This posistion is not devined in the list
            #so for practicality place the most outer 3 (3 axis) coordinates in the last list
            if pos1 == windowx:
                pos1 = windowx-1
            elif pos2 == windowy:
                pos2 = windowy-1
            elif pos3 ==windowz:
                pos3 = windowz-1
            
            binxyz[pos1][pos2][pos3].append(data_order[i])
            number_processed +=1

print("number of designated particles",number_processed)
count = 0
maximum = 0 
averaging_field = []
gauss_field = []
for i in range(0,len(binxyz)):
    for k in range(0,len(binxyz[i])):
        for m in range(0,len(binxyz[i][k])):
            #check howmany bins are empty
            if binxyz[i][k][m]==[]:
                    count +=1
            #for every bin sum the velocities and devide by the amount of particles inside the bin
            elif binxyz[i][k][m]!=[]:
                u_count,v_count,w_count = 0,0,0
                u_gsum,v_gsum,w_gsum = 0,0,0
                for n in range(0,len(binxyz[i][k][m])):
                    u_count = u_count + binxyz[i][k][m][n][3]
                    v_count = u_count + binxyz[i][k][m][n][4]
                    w_count = u_count + binxyz[i][k][m][n][5]
                u_ave,v_ave,w_ave = u_count/len(binxyz[i][k][m]),v_count/len(binxyz[i][k][m]),w_count/len(binxyz[i][k][m])

                #This is all for the comments part
                for n in range(0, len(binxyz[i][k][m])):
                    u_gsum = u_gsum + ((binxyz[i][k][m][n][3] - u_ave) ** 2)
                    v_gsum = v_gsum + ((binxyz[i][k][m][n][4] - v_ave) ** 2)
                    w_gsum = w_gsum + ((binxyz[i][k][m][n][5] - w_ave) ** 2)
                u_sd,v_sd,w_sd=math.sqrt(u_gsum / len(binxyz[i][k][m])), math.sqrt(v_gsum / len(binxyz[i][k][m])), math.sqrt(w_gsum / len(binxyz[i][k][m]))
                for n in range(0, len(binxyz[i][k][m])):
                    u_weight = (1/(u_sd*math.sqrt(2*math.pi)))*math.exp((-1/2)*(((binxyz[i][k][m][3]-u_ave)/u_sd)**2))*binxyz[i][k][m][3]
                    v_weight = (1/(v_sd*math.sqrt(2*math.pi)))*math.exp((-1/2)*(((binxyz[i][k][m][4]-v_ave)/v_sd)**2))*binxyz[i][k][m][4]
                    w_weight = (1/(w_sd*math.sqrt(2*math.pi)))*math.exp((-1/2)*(((binxyz[i][k][m][5]-w_ave)/u_sd)**2))*binxyz[i][k][m][5]

                #use the location list to establish the mid point of the bin
                averaging_field.append([(x_loc[i]+x_loc[i+1])/2,(y_loc[k]+y_loc[k+1])/2,(z_loc[m]+z_loc[m+1])/2,u_ave,v_ave,w_ave])
                gauss_field.append([(x_loc[i]+x_loc[i+1])/2,(y_loc[k]+y_loc[k+1])/2,(z_loc[m]+z_loc[m+1])/2,u_gauss,v_gauss,w_gauss])
            #check for the bin with the maximum amount of particles
            if len(binxyz[i][k][m])>maximum:
                maximum=len(binxyz[i][k][m])
print("Number of bins",windowx*windowy*windowz)
print("Number of empty bins",count)
print("Maximum amount of particles in bin",maximum)

print(gauss_field[0:4])


#define the dimensions of the measurmement volume
print("x:", x_min,x_max)
print("y:", y_min,y_max)
print("z:", z_min,z_max)

print("Volume:",(x_max-x_min)*(y_max-y_min)*(z_max-z_min)*10**(-9),"m^3")
print("Sub-volume",delta_x*delta_y*delta_z,"mm^3")

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#nrParticles = 500
#dataSel = data[0:nrParticles,:]
#ax.quiver(dataSel[:,0],dataSel[:,1],dataSel[:,2],dataSel[:,3],dataSel[:,4],dataSel[:,5])
#plt.show()
