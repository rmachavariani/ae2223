import numpy as np
import time
import pdb
import statistics as sts



def read_cell():
    ''' Loads data from carMirrorData.dat
    :return: nrVectors x 6 np array
    '''
    # load the data
    t1 = time.time()
    data = np.loadtxt("example.dat",max_rows = None, delimiter=",")
    t2 = time.time()
    print("Loading done in ", "{:.2f}".format(t2 - t1), " s")

    return data





def detect_outliers(particle_data):
    data = []
    st_dev = sts.stdev(particle_data[:,3])
    for i in range(np.size(particle_data, axis = 0)):
        if abs(particle_data[i,3] - np.mean(particle_data[:,3])) < st_dev:
            data.append(particle_data[i,0:6])

    return np.array(data)


inp = read_cell()
a = detect_outliers(inp)
print(a)
#pdb.set_trace()

