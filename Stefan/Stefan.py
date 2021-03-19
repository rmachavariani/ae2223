import timing
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np
import threading
import time


#####################################################################
#                                                                   #
#                       Standard functions                          #
#                                                                   #
#####################################################################

def upload():
    t1 = time.time()
    # Define the structure of the data
    dtype1 = np.dtype([("x", float), ("y", float), ("z", float), ("u", float), ("v", float), ("w", float)])
    datalocation = "/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/My coding/carMirrorData.dat"
    # upload the data from the file to a numpy array with designated structure
    data = np.loadtxt(datalocation, dtype=dtype1)

    # order the data values with x from small to large x
    data_order = np.sort(data, order="x")
    t2 = time.time()
    print("Number of particles", len(data_order), "orderd in:", round(t2 - t1, 3), "sec")
    return data, data_order, dtype1


def bincheck(binxyz):
    count = 0
    for i in range(0, len(binxyz)):
        for k in range(0, len(binxyz[i])):
            for m in range(0, len(binxyz[i][k])):
                # check how many bins are empty
                if binxyz[i][k][m] == []:
                    count += 1
    return count


def histogram(binxyz, binx, biny, binz, dtype1):
    lst = np.array(binxyz[binx][biny][binz], dtype=dtype1)
    order = ["u", "v", "w"]

    plt.figure()
    d_v = np.linspace(min(min(lst["u"]), min(lst["v"]), min(lst["w"])) + 10,
                      max(max(lst["u"]), max(lst["v"]), max(lst["w"])) - 10, 30, dtype=np.float64)

    plt.subplot(1, 1, 1)
    plt.title("Histogram bin:" + str(binx) + "," + str(biny) + "," + str(binz))
    plt.hist([lst["u"], lst["v"], lst["w"]], bins=d_v, label=["u", "v", "w"])
    plt.legend(loc='upper right')
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Count")
    plt.show()


#####################################################################
#                                                                   #
#                            Classes                                #
#                                                                   #
#####################################################################

class GridBin:
    def __init__(self, windowx, windowy, windowz, offset, data, data_order):
        self.offset = offset
        self.windowx = windowx
        self.windowy = windowy
        self.windowz = windowz
        self.data = data
        self.data_order = data_order

    def boundaries(self):
        # define outer bounds of the measurement volume
        x_min, x_max = min(self.data["x"]), max(self.data["x"])
        y_min, y_max = min(self.data["y"]), max(self.data["y"])
        z_min, z_max = min(self.data["z"]), max(self.data["z"])

        return x_min, x_max, y_min, y_max, z_min, z_max

    def bins(self, offset, x_min, x_max, y_min, y_max, z_min, z_max):
        # can make a offset were to start compared to the min and maximum location at the axis
        delta_x = ((x_max + self.offset[0][1]) - (x_min + self.offset[0][0])) / self.windowx
        delta_y = ((y_max + self.offset[1][1]) - (y_min + self.offset[1][0])) / self.windowy
        delta_z = ((z_max + self.offset[2][1]) - (z_min + self.offset[2][0])) / self.windowz

        # bound locations of all the bins (use of +1 because of np.digitize in distribute function)
        x_bound = np.linspace(x_min + self.offset[0][0], x_max + self.offset[0][1], windowx + 1, dtype=np.float64)
        y_bound = np.linspace(y_min + self.offset[1][0], y_max + self.offset[1][1], windowy + 1, dtype=np.float64)
        z_bound = np.linspace(z_min + self.offset[2][0], z_max + self.offset[2][1], windowz + 1, dtype=np.float64)

        # List that define the center locations of every bin
        x_loc = np.arange(start=x_min + self.offset[0][0] + 0.5 * delta_x, stop=x_max + self.offset[0][1], step=delta_x)
        y_loc = np.arange(start=y_min + self.offset[1][0] + 0.5 * delta_y, stop=y_max + self.offset[1][1], step=delta_y)
        z_loc = np.arange(start=z_min + self.offset[2][0] + 0.5 * delta_z, stop=z_max + self.offset[2][1], step=delta_z)

        return delta_x, delta_y, delta_z, x_bound, y_bound, z_bound, x_loc, y_loc, z_loc

    def grid(self):
        # make list with all bins were the particle information can be added to
        binxyz = [[[[] for k in range(self.windowz)] for l in range(self.windowy)] for i in
                  range(self.windowx)]  # [x:[y:[z:[particle data:()]]]]

        return binxyz

    def distribute(self, binxyz, x_bound, y_bound, z_bound):
        t1 = time.time()
        number_processed = 0

        # np.digitize compares the location of the particle with the boundaries of the bins
        # and returns a list with the index location of the lowest boundarie of the bin
        x_in = np.digitize(data_order["x"], x_bound)
        y_in = np.digitize(data_order["y"], y_bound)
        z_in = np.digitize(data_order["z"], z_bound)

        # Placing all the particles in the correct place in binxyz
        for i in range(len(data_order)):
            # When a particles coordinate is smaller then the smallest boundarie it is given 0. When
            # larger then the largest boundarie it is given len(bound). So using if statement
            # particles outside the designated measurement volume are excluded
            if 1 <= x_in[i] < len(x_bound) and 1 <= y_in[i] < len(y_bound) and 1 <= z_in[i] < len(z_bound):
                # use -1 as indecis are starting at 1 and list start at 0
                binxyz[x_in[i] - 1][y_in[i] - 1][z_in[i] - 1].append(data_order[i])
                # count the number of processed particles
                number_processed += 1

        t2 = time.time()
        print("Number of designated particles", number_processed, "in:", round(t2 - t1, 3), "sec")
        return binxyz, number_processed


class Averaging:
    def __init__(self, binxyz, x_loc, y_loc, z_loc, dtype1):
        self.binxyz = binxyz
        self.x_loc = x_loc
        self.y_loc = y_loc
        self.z_loc = z_loc
        self.dtype1 = dtype1

    def mean(self):
        t1 = time.time()
        maximum = 0
        averaging_field = []
        for i in range(0, len(self.binxyz)):
            for k in range(0, len(self.binxyz[i])):
                for m in range(0, len(self.binxyz[i][k])):
                    # for every bin sum the velocities and devide by the amount of particles inside the bin
                    # can only be done for bins with particles inside
                    if self.binxyz[i][k][m] != []:
                        # Convert the particle information in the bin to an array for more simple use
                        lst = np.array(self.binxyz[i][k][m], dtype=self.dtype1)
                        # Use numpy for the mean value so that no loop has to be used
                        u_ave = np.average(lst["u"])
                        v_ave = np.average(lst["v"])
                        w_ave = np.average(lst["w"])
                        # append the mean average velocity components to a list with center location of bin
                        averaging_field.append([self.x_loc[i], self.y_loc[k], self.z_loc[m], u_ave, v_ave, w_ave])
                    # check for the bin with the maximum amount of particles
                    if len(self.binxyz[i][k][m]) > maximum:
                        maximum = len(self.binxyz[i][k][m])

        t2 = time.time()
        print("Mean averaging takes", round(t2 - t1, 3), "sec")

        return maximum, averaging_field, u_ave, v_ave, w_ave

    def gaussian(self):
        t1 = time.time()
        gauss_field = []
        for i in range(0, len(self.binxyz)):
            for k in range(0, len(self.binxyz[i])):
                for m in range(0, len(self.binxyz[i][k])):
                    # for every bin sum the velocities and devide by the amount of particles inside the bin
                    if self.binxyz[i][k][m] != []:
                        u_gauss, v_gauss, w_gauss = 0, 0, 0
                        u_weight, v_weight, w_weight = 0, 0, 0
                        u_sum, v_sum, w_sum = 0, 0, 0

                        # If the number of particles inside a bin is equal to 1 then sd=0 and
                        # calculating the normal distrubtion results in a devision by 0
                        if len(binxyz[i][k][m]) > 1:
                            # make the list a array so that information can be extracted easily
                            lst = np.array(self.binxyz[i][k][m], dtype=self.dtype1)
                            # use numpy to calculate the mean and standard deviation
                            u_ave, u_sd = np.average(lst["u"]), np.std(lst["u"])
                            v_ave, v_sd = np.average(lst["v"]), np.std(lst["v"])
                            w_ave, w_sd = np.average(lst["w"]), np.std(lst["w"])

                            for n in range(0, len(lst)):
                                # calculate for every velocity component the normal value and sum this together for the bin
                                u_nor = (1 / (u_sd * np.sqrt(2 * np.pi))) * np.exp(
                                    (-1 / 2) * (((lst[n][3] - u_ave) / u_sd) ** 2))
                                u_weight = u_weight + u_nor * lst[n][3]
                                u_sum = u_sum + u_nor

                                v_nor = (1 / (v_sd * np.sqrt(2 * np.pi))) * np.exp(
                                    (-1 / 2) * (((lst[n][4] - v_ave) / v_sd) ** 2))
                                v_weight = v_weight + v_nor * lst[n][4]
                                v_sum = v_sum + v_nor

                                w_nor = (1 / (w_sd * np.sqrt(2 * np.pi))) * np.exp(
                                    (-1 / 2) * (((lst[n][5] - w_ave) / w_sd) ** 2))
                                w_weight = w_weight + w_nor * lst[n][5]
                                w_sum = w_sum + w_nor

                            # devide the sum of the normal*velocity by the sum of the normal values to get the gaussian velocity
                            u_gauss = u_weight / (u_sum)
                            v_gauss = v_weight / (v_sum)
                            w_gauss = w_weight / (w_sum)

                        # when there is only 1 particle in a bin this is automaticly the determining particle velocity
                        elif len(binxyz[i][k][m]) == 1:
                            u_gauss = self.binxyz[i][k][m][0][3]
                            v_gauss = self.binxyz[i][k][m][0][4]
                            w_gauss = self.binxyz[i][k][m][0][5]

                        # use the location list to establish the mid point of the bin
                        gauss_field.append([self.x_loc[i], self.y_loc[k], self.z_loc[m], u_gauss, v_gauss, w_gauss])
        t2 = time.time()
        print("Gaussian averaging takes", round(t2 - t1, 3), "sec")
        return gauss_field

    def polynomialfit(self):
        polynomial_field = []
        t1 = time.time()
        for i in range(0, len(self.binxyz)):
            for k in range(0, len(self.binxyz[i])):
                for m in range(0, len(self.binxyz[i][k])):
                    if len(binxyz[i][k][m]) >= 3:
                        lst = np.array(self.binxyz[i][k][m], dtype=self.dtype1)

                        for axis in range(0, 3):
                            if axis ==0:
                                points = lst["x"]
                                velocity = lst["u"]
                            elif axis ==1:
                                points = lst["y"]
                                velocity = lst["v"]
                            elif axis ==2:
                                points = lst["z"]
                                velocity = lst["w"]

                            A = np.vstack([np.ones(len(points)),points,(points**2)]).T
                            s,d,f = np.linalg.lstsq(A,velocity, rcond=None)[0]

                            if axis ==0:
                                u_pol = s+d*self.x_loc[i]+f*(self.x_loc[i]**2)
                                if u_pol >= 15:
                                    print("u",u_pol,i,k,m)
                                    print(s,d,f)
                                    print(lst["u"])
                            elif axis ==1:
                                v_pol = s+d*self.y_loc[k]+f*(self.y_loc[k]**2)
                                if v_pol >= 5:
                                    print("v",v_pol,i,k,m)
                                    print(s,d,f)
                                    print(lst["v"])
                            elif axis ==2:
                                w_pol = s+d*self.z_loc[m]+f*(self.z_loc[m]**2)
                                if w_pol >= 5:
                                    print("w",w_pol,i,k,m)
                                    print(s,d,f)
                                    print(lst["w"])

                        polynomial_field.append([x_loc[i], y_loc[k], z_loc[m], u_pol, v_pol, w_pol])
        t2 = time.time()
        print("Polynomial fit done in", round(t2 - t1, 3), "sec")
        return polynomial_field


#####################################################################
#                                                                   #
#                           Main script                             #
#                                                                   #
#####################################################################

##########################Import data################################

data, data_order, dtype1 = upload()

#######################Setting up the grid###########################

# define measurement volume parameters
windowx, windowy, windowz = 10, 10, 10
offset = [[0, 0], [0, 0], [0, 0]]
particle = GridBin(windowx, windowy, windowz, offset, data, data_order)

# calculate the maximum and minimum coordinates in all 3 directions
x_min, x_max, y_min, y_max, z_min, z_max = particle.boundaries()

# amount of bins that are placed along the axis
delta_x, delta_y, delta_z, x_bound, y_bound, z_bound, x_loc, y_loc, z_loc = particle.bins(offset, x_min, x_max, y_min,
                                                                                          y_max, z_min, z_max)

# create the desired grid in which the particles have to be devided
binxyz = particle.grid()

#######################Filling in the grid###########################

# apply all the particles that are located within the designated measurement volume to a bin
binxyz, number_processed = particle.distribute(binxyz, x_bound, y_bound, z_bound)

count = bincheck(binxyz)
print("Number of bins", windowx * windowy * windowz)
print("Number of empty bins", count, "(", round(count / (windowx * windowy * windowz), 2) * 100, "%)")

############################Averaging#################################

general = Averaging(binxyz, x_loc, y_loc, z_loc, dtype1)

maximum, averaging_field, u_ave, v_ave, w_ave = general.mean()

gauss_field = general.gaussian()

polynomial_field = general.polynomialfit()

print("Maximum amount of particles in bin", maximum)

# define the dimensions of the measurmement volume
print("x:", x_min, x_max)
print("y:", y_min, y_max)
print("z:", z_min, z_max)

print("Volume:", (x_max - x_min) * (y_max - y_min) * (z_max - z_min) * 10 ** (-9), "m^3")
print("Sub-volume", delta_x * delta_y * delta_z, "mm^3")

# plot the histogram of a certain bin
binx, biny, binz = int(windowx / 2), int(windowy / 2), int(windowz / 2)
histogram(binxyz, binx, biny, binz, dtype1)

#####################################################################
#                                                                   #
#                           Plot script                             #
#                                                                   #
#####################################################################

fig = plt.figure()
###################
# Normal plot
###################
ax = fig.add_subplot(1, 3, 1, projection='3d')
plt.title("Averaging")

x, y, z, u, v, w = zip(*averaging_field)

average_plt = ax.quiver(x, y, z, u, v, w, length=3)  # , color="blue", colormap(norm(colors)), cmap=cm.coolwarm
# fig.colorbar(average_plt, orientation='vertical')

ax.set_xlim([x_min, x_max])
ax.set_ylim([y_min, y_max])
ax.set_zlim([z_min, z_max])
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('z [mm]')

###################
# Guassian plot
###################
ax = fig.add_subplot(1, 3, 2, projection='3d')
plt.title("Gaussian")

x, y, z, u, v, w = zip(*gauss_field)

gaussian_plt = ax.quiver(x, y, z, u, v, w, length=3)  # , color="blue", colormap(norm(colors)), cmap=cm.coolwarm
# fig.colorbar(gaussian_plt, orientation='vertical')

ax.set_xlim([x_min, x_max])
ax.set_ylim([y_min, y_max])
ax.set_zlim([z_min, z_max])
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('z [mm]')

###################
# Polynomial plot
###################

ax = fig.add_subplot(1, 3, 3, projection='3d')
plt.title("Polynomial")

x, y, z, u, v, w = zip(*polynomial_field)

polynomial_plt = ax.quiver(x, y, z, u, v, w, length=3)  # , color="blue", colormap(norm(colors)), cmap=cm.coolwarm
# fig.colorbar(gaussian_plt, orientation='vertical')

ax.set_xlim([x_min, x_max])
ax.set_ylim([y_min, y_max])
ax.set_zlim([z_min, z_max])
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('z [mm]')

###################
# Bin visualisation
###################
# ax = fig.add_subplot(1, 3, 3, projection='3d')
# plt.title("Bin")
# test =binxyz[int(windowx/2)][int(windowy/2)][int(windowz/2)]
# x,y,z,u,v,w = zip(*test)

# # colors = np.arctan2(u, v)
# # norm = Normalize()
# # norm.autoscale(colors)
# # colormap = cm.coolwarm

# gaussian_plt=ax.quiver(x, y, z, u, v, w, length=0.1) #c, color="blue", olormap(norm(colors)), cmap=cm.coolwarm
# #fig.colorbar(gaussian_plt, orientation='vertical')

# ax.set_xlim([x_bound[int(windowx/2)], x_bound[int(windowx/2)+1]])
# ax.set_ylim([y_bound[int(windowy/2)], y_bound[int(windowy/2)+1]])
# ax.set_zlim([z_bound[int(windowz/2)], z_bound[int(windowz/2)+1]])
# ax.set_xlabel('x [mm]')
# ax.set_ylabel('y [mm]')
# ax.set_zlabel('z [mm]')

plt.show()