from os import listdir
from PIL import Image, ImageFilter
import numpy as np
from matplotlib import pyplot
from preprocessing import *




#testimg = None
#counter = 0
#take = 3

# Dimension of all images
imgx = 300
imgy = 400

def getInputData(imgx, imgy):

    pictures = []

    # Fetch the pictures
    for picName in listdir("data"):
        picture = Image.open("data/" + picName)
        pictures.append(picture)

    # Initialize the output list
    input_list = []
    for image in pictures:

        #print(image.size)

        image = image.convert("L")                      # Convert to grayscale
        image = image.filter(ImageFilter.FIND_EDGES)    # Apply edgedetect filter

        #------ code for testing -------
        #if counter == take:
        #    testimg = image
        #pix = np.array(image)
        #pix = np.array(image.getdata()).reshape(image.size[0], image.size[1], 3)


        pix = np.array(image.resize((imgx, imgy)))      # Resize to set dimension
        pix = pix.reshape(-1,1)                         # Flatten
        input_list.append(pix)

        #------ code for testing -------
        #counter += 1
    return np.array(input_list)

#---- code for testing ------
stuff = getInputData(imgx, imgy)
print(stuff[1])

#pyplot.imshow(testimg)
#pyplot.show()

#print(input_list[1][10,10])



