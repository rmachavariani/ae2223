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

        # Flatten the RGB image
        image_array = np.array(image.resize((imgx, imgy)))                 # Convert image to array
        image_flat = image_array.reshape(-1,1)                             # Flatten

        # Prepare data using edgedetect filter in greyscale
        image_grey = image.convert("L")                                    # Convert to grayscale
        image_grey_filtered = image_grey.filter(ImageFilter.FIND_EDGES)    # Apply edgedetect filter
        pix = np.array(image_grey_filtered.resize((imgx, imgy)))           # Resize to set dimension
        data_gf = pix.reshape(-1,1)                                        # Flatten

        pix = data_gf + image_flat                                         # Concatenate the grey and rgb image arrays

        # Check if class of type exists, if so add data to class, if not create class of type
        if picName in names:
                picName.data.append(pix)
        else:
            input_list.append(pix)
            names.append(picName)
            i = names.index(picName)
            input_list[i] = testObject(pix, picName)

        return input_list

#---- code for testing ------
#stuff = getInputData(imgx, imgy)
#print(stuff[0].actype)





