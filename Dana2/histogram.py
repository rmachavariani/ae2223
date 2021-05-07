import matplotlib.pyplot as plt
from Dana2.class_def import *
import seaborn as sns
from Dana2.gridConstructionSphereFast import getSphericalGridWithVectorsFast



def getMostPopulated(data, n_bins):
    print('getting the most populated')
    sizeX = np.size(data,axis=0)
    sizeY = np.size(data,axis=1)
    sizeZ = np.size(data,axis=2)

    lenArray = []

    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                lenArray.append([i, j, k, len(data[i][j][k].vectors)])


    lenArray = np.array(lenArray)
    largestKeys = np.argpartition(lenArray[:,3], -n_bins)[-n_bins:]
    bins = lenArray[largestKeys]

    return bins



def getHistograms(data, bins):

    for bin in bins:
        bin = data[bin[0], bin[1], bin[2]]
        vectors = bin.vectors
        u_values = []
        v_values = []
        w_values = []

        for vector in vectors:
            u_values.append(vector.u)
            v_values.append(vector.v)
            w_values.append(vector.w)

        sns.set_theme()
        sns.distplot(u_values, hist=False, kde=True,
                     kde_kws={'linewidth': 1},
                     label='u')
        sns.distplot(v_values, hist=False, kde=True,
                     kde_kws={'linewidth': 1},
                     label='v')
        sns.distplot(w_values, hist=False, kde=True,
                     kde_kws={'linewidth': 1},
                     label='w')

        # plt.hist(u_values, bins='auto', alpha=0.5, label="U", color="red", rwidth=0.85, histtype='step', fc='none', lw=1.5)
        # plt.hist(v_values,  bins='auto', alpha=0.5, label="V",  color="blue", rwidth=0.85, histtype='step', fc='none', lw=1.5)
        # plt.hist(w_values,  bins='auto', alpha=0.5, label="W", color="lime",  rwidth=0.85, histtype='step', fc='none', lw=1.5)
        plt.legend(loc='upper left')

        plt.xlabel("Velocity [m/s]")
        plt.show()



n_bins = 4
# data = getRectangularGridWithVectors(5000, 5000, 5000)
data = getSphericalGridWithVectorsFast(50,50, None)
bins = getMostPopulated(data, n_bins)
# getGeneralDensity()
# histograms2d = get2dHistograms(data,bins)
histograms = getHistograms(data, bins)