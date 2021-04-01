from os import listdir
from PIL import Image, ImageFilter
import numpy as np
from matplotlib import pyplot
from preprocessing import *

# Dimension of all images
imgx = 300
imgy = 400


class testObject:
    def __init__(self, data, type):
        self.data = [data]
        self.actype = type


def getInputData(imgx, imgy):
    '''
    Call function to get the preprocessed data

    Returns
    -------
    input_list : List of testObjects
        a testObject contains:
         - 'data' : list which contains lists of flattened, concatenated (edgedetected) grey and RGB image data.
        One for each picture that is in the dataset.
         - 'actype' : string of the aircraft type.
    '''
    # Initiale list for output
    input_list = []
    names = []

    for picName in listdir("data"):
        image = Image.open("data/" + picName)

        # Prepare name
        picName = picName.split(".")[0]                                    # Delete file extension
        picNameSplit = picName.split("_")                                  # Split filename at underscore
        picNameSplit.remove(picNameSplit[-1])                              # Remove the picture number
        picName = "-".join(picNameSplit)                                   # Join name together with hyphen
        print(picName)
        # Flatten the RGB image
        image_array = np.array(image.resize((imgx, imgy)))                 # Convert image to array
        image_flat = image_array.reshape(-1,1)                             # Flatten

        # Prepare data using edgedetect filter in greyscale
        image_grey = image.convert("L")                                    # Convert to grayscale
        image_grey_filtered = image_grey.filter(ImageFilter.FIND_EDGES)    # Apply edgedetect filter
        pix = np.array(image_grey_filtered.resize((imgx, imgy)))           # Resize to set dimension
        data_gf = pix.reshape(-1,1)                                        # Flatten

        pix = np.concatenate((data_gf, image_flat))                                         # Concatenate the grey and rgb image arrays

        # Check if class of type exists, if so add data to class, if not create class of type
        if picName in names:
            i = names.index(picName)
            input_list[i].data.append(pix)
        else:
            input_list.append(pix)
            names.append(picName)
            i = names.index(picName)
            input_list[i] = testObject(pix, picName)

    return input_list, names




def create_arrays(imgx, imgy):
    data_list, names = getInputData(imgx, imgy)
    typeqty = len(data_list)
    dataqty = 0
    imagesize = len(data_list[0].data[0])
    for ob in data_list:
        datasize = len(ob.data)
        dataqty += datasize

    data_array = np.zeros((imagesize,dataqty))
    type_array = np.zeros((typeqty, dataqty))

    i = 0
    for typedata in data_list:
        for imagedata in typedata.data:
            data_array[:,i] = np.reshape(imagedata,imagesize)
            j = names.index(typedata.actype)
            type_array[j,i] = 1
            i += 1
    return data_array, type_array, names

#---- code for testing ------
#input_list, names = getInputData(imgx, imgy)
#data_array, type_array, names = create_arrays(input_list, names)

#print(data_array.shape, type_array, names)
