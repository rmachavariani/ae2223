import timing
import matplotlib.pyplot as plt
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
    # Import average monthly precip to numpy array
    dtype1 = np.dtype([("x", float), ("y", float), ("z", float), ("u", float), ("v", float), ("w", float)])
    # datalocation = "/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/Coding/carMirrorData.dat"
    datalocation = "/Users/stefanrooze/Documents/TU Delft/Quarter 3/AE2223-I Test analysees & Simulation/My coding/carMirrorData.dat"
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
                # check howmany bins are empty
                if binxyz[i][k][m] == []:
                    count += 1
    return count


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

        # List that define the center location locations of every bin
        x_loc = np.arange(start=x_min + self.offset[0][0] + 0.5 * delta_x, stop=x_max + self.offset[0][1], step=delta_x)
        y_loc = np.arange(start=y_min + self.offset[1][0] + 0.5 * delta_y, stop=y_max + self.offset[1][1], step=delta_y)
        z_loc = np.arange(start=z_min + self.offset[2][0] + 0.5 * delta_z, stop=z_max + self.offset[2][1], step=delta_z)

        return delta_x, delta_y, delta_z, x_loc, y_loc, z_loc

    def grid(self):
        # make list with all bins were the particle information can be added to
        binxyz = [[[[] for k in range(self.windowz)] for l in range(self.windowy)] for i in
                  range(self.windowx)]  # [x:[y:[z:[particle data:()]]]]

        return binxyz

    def distribute(self, binxyz, x_min, x_max, y_min, y_max, z_min, z_max, delta_x, delta_y, delta_z):
        # with the use of dividing the location of the particles by the bin width
        # the rounded number results in the bin the particle has to be placed in
        t1 = time.time()
        number_processed = 0
        for i in range(len(data_order)):
            # make sure that the particle is inside the measurement volume defined with offset
            if self.data_order[i][0] + abs(x_min) >= self.offset[0][0] and self.data_order[i][1] + abs(y_min) >= \
                    self.offset[1][0] and self.data_order[i][2] + abs(z_min) >= self.offset[2][0]:
                if self.data_order[i][0] <= x_max + self.offset[0][1] and self.data_order[i][1] <= y_max + \
                        self.offset[1][1] and self.data_order[i][2] <= z_max + self.offset[2][1]:
                    # as there are negative coordinate points add the minimum value to get positive values
                    # that can be used to determine location
                    x = self.data_order[i][0] - x_min - self.offset[0][0]
                    # int makes sure that the particle is assigned to the correct place in the list
                    pos1 = int(x / delta_x)

                    y = self.data_order[i][1] - y_min - self.offset[1][0]
                    pos2 = int(y / delta_y)

                    z = self.data_order[i][2] - z_min - self.offset[2][0]
                    pos3 = int(z / delta_z)

                    # the maximum coordinates will be assigned pos=window. This posistion is not devined in the list
                    # so for practicality place the most outer 3 (3 axis) coordinates in the last list
                    if pos1 == self.windowx:
                        pos1 = self.windowx - 1
                    if pos2 == self.windowy:
                        pos2 = self.windowy - 1
                    if pos3 == self.windowz:
                        pos3 = self.windowz - 1

                    binxyz[pos1][pos2][pos3].append(self.data_order[i])
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
                        lst = np.array(self.binxyz[i][k][m], dtype=self.dtype1)
                        u_ave = np.mean(lst["u"])
                        v_ave = np.mean(lst["v"])
                        w_ave = np.mean(lst["w"])
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
                        u_lst, v_lst, w_lst = [], [], []
                        u_gauss, v_gauss, w_gauss = 0, 0, 0
                        u_weight, v_weight, w_weight = 0, 0, 0

                        # If the number of particles inside a bin is equal to 1 then sd=0 and
                        # calculating the normal distrubtion results in a devision by 0
                        if len(binxyz[i][k][m]) > 1:
                            # make the list a array so that information can be extracted easily
                            lst = np.array(self.binxyz[i][k][m], dtype=self.dtype1)
                            # use numpy to calculate the mean and standard deviation
                            u_ave, u_sd = np.mean(lst["u"]), np.std(lst["u"])
                            v_ave, v_sd = np.mean(lst["v"]), np.std(lst["v"])
                            w_ave, w_sd = np.mean(lst["w"]), np.std(lst["w"])

                            for n in range(0, len(lst)):
                                u_nor = (1 / (u_sd * np.sqrt(2 * np.pi))) * np.exp(
                                    (-1 / 2) * (((lst[n][3] - u_ave) / u_sd) ** 2))
                                u_weight = u_weight + u_nor * lst[n][3]
                                u_lst.append(u_nor)

                                v_nor = (1 / (v_sd * np.sqrt(2 * np.pi))) * np.exp(
                                    (-1 / 2) * (((lst[n][4] - v_ave) / v_sd) ** 2))
                                v_weight = v_weight + v_nor * lst[n][4]
                                v_lst.append(v_nor)

                                w_nor = (1 / (w_sd * np.sqrt(2 * np.pi))) * np.exp(
                                    (-1 / 2) * (((lst[n][5] - w_ave) / w_sd) ** 2))
                                w_weight = w_weight + w_nor * lst[n][5]
                                w_lst.append(w_nor)

                            u_gauss = u_weight / (sum(u_lst))
                            v_gauss = v_weight / (sum(v_lst))
                            w_gauss = v_weight / (sum(v_lst))

                        elif len(binxyz[i][k][m]) == 1:
                            u_gauss = self.binxyz[i][k][m][0][3]
                            v_gauss = self.binxyz[i][k][m][0][4]
                            w_gauss = self.binxyz[i][k][m][0][5]

                        # use the location list to establish the mid point of the bin
                        gauss_field.append([self.x_loc[i], self.y_loc[k], self.z_loc[m], u_gauss, v_gauss, w_gauss])
        t2 = time.time()
        print("Gaussian averaging takes", round(t2 - t1, 3), "sec")
        return gauss_field


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
offset = [[100, -100], [100, -100], [100, -100]]
particle = GridBin(windowx, windowy, windowz, offset, data, data_order)

# calculate the maximum and minimum coordinates in all 3 directions
x_min, x_max, y_min, y_max, z_min, z_max = particle.boundaries()

# amount of bins that are placed along the axis
delta_x, delta_y, delta_z, x_loc, y_loc, z_loc = particle.bins(offset, x_min, x_max, y_min, y_max, z_min, z_max)

# create the desired grid in which the particles have to be devided
binxyz = particle.grid()

#######################Filling in the grid###########################

# apply all the particles that are located within the designated measurement volume to a bin
binxyz, number_processed = particle.distribute(binxyz, x_min, x_max, y_min, y_max, z_min, z_max, delta_x, delta_y,
                                               delta_z)

count = bincheck(binxyz)
print("Number of bins", windowx * windowy * windowz)
print("Number of empty bins", count)

############################Averaging#################################

general = Averaging(binxyz, x_loc, y_loc, z_loc, dtype1)

maximum, averaging_field, u_ave, v_ave, w_ave = general.mean()

gauss_field = general.gaussian()
print("Maximum amount of particles in bin", maximum)

# define the dimensions of the measurmement volume
print("x:", x_min, x_max)
print("y:", y_min, y_max)
print("z:", z_min, z_max)

print("Volume:", (x_max - x_min) * (y_max - y_min) * (z_max - z_min) * 10 ** (-9), "m^3")
print("Sub-volume", delta_x * delta_y * delta_z, "mm^3")

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# nrParticles = 500
# dataSel = data[0:nrParticles,:]
# ax.quiver(dataSel[:,0],dataSel[:,1],dataSel[:,2],dataSel[:,3],dataSel[:,4],dataSel[:,5])
# plt.show()


# x_min,x_max,y_min,y_max,z_min,z_max = boundaries(data)
# #amount of bins that are placed along the axis
# windowx, windowy, windowz = 10, 10, 10
# offset = [[100,-100],[100,-100],[100,-100]]
# delta_x,delta_y,delta_z,x_loc,y_loc,z_loc = bins(data_order,windowx,windowy,windowz,offset,x_min,x_max,y_min,y_max,z_min,z_max)
# binxyz = grid(windowx,windowy,windowz)