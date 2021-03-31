from os import listdir
from PIL import Image, ImageFilter
import numpy as np
from matplotlib import pyplot
from preprocessing import *

# Dimension of all images
imgx = 300
imgy = 400

def getInputData(imgx, imgy):

    # Initiale list for output
    input_list = []

    for picName in listdir("data"):
        image = Image.open("data/" + picName)
        image = image.convert("L")                      # Convert to grayscale
        image = image.filter(ImageFilter.FIND_EDGES)    # Apply edgedetect filter
        pix = np.array(image.resize((imgx, imgy)))      # Resize to set dimension
        pix = pix.reshape(-1,1)                         # Flatten
        input_list.append(pix)

    return np.array(input_list)

#---- code for testing ------
stuff = getInputData(imgx, imgy)
print(stuff[1])





