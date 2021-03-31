from os import listdir
from PIL import Image
import numpy as np

def loadRawData():
    
    # empty list
    pictures = list()
    names = list()
    
    # loop over pictures in data
    for picName in listdir("data"):
        
        # load every picture as a numpy array and add to list
        picture = Image.open("data/"+picName)
        pictures.append(picture)
        
        # process pictures name and add to list
        picName = picName.split(".")[0]
        picNameSplit = picName.split("_")
        picNameSplit.remove(picNameSplit[-1])
        picName = "-".join(picNameSplit)
        
        names.append(picName)
    
    return pictures, names


def resizePictures(pictures):
    
    # find the smallest dimensions
    smallx = 10000
    smally = 10000
    
    # look for the smallest dimension in each direction
    for picture in pictures:
        
        shapex, shapey = picture.size
        
        if shapex < smallx:
            
            smallx = shapex
        
        if shapey < smally:
            
            smally = shapey
    
    # reshape all pictures
    for i in range(len(pictures)):
        
        # resizing all pictures to same dimensions
        pictures[i] = np.array(pictures[i].resize((smallx,smally)))
    
    
    return pictures
        
def flattenPictures(pictures):
    
    for i in range(len(pictures)):
        
        pictures[i] = pictures[i].reshape(-1,1)
    
    return pictures


#------------------------------------------------------------#
#------------------------------------------------------------#
#------------------------------------------------------------#
#----------------------------MAIN----------------------------#
#------------------------------------------------------------#
#------------------------------------------------------------#
#------------------------------------------------------------#



def getPreprocessedData():
    '''
    Call function to get the preprocessed data

    Returns
    -------
    labels : List of strings
        Type of airplane belonging to each picture.
        
    pictures : list of pictures as used in PIL
        Original pictures in same order as labels.
        
    shaped_pictures : List of numpy arrays
        Reshaped pictures (all same dimension) 
        and transformed into a numpy array.
        
    flat_pictures : List of numpy arrays
        Pictures are flattened into numpy arrays with shape (nrPixels x 1)

    '''
    
    # pictures is list of images, names is list of strings
    pictures, labels = loadRawData()
    picturesBis = pictures

    # pictures is now converted to list of numpy arrays with same dimensions
    shaped_pictures = resizePictures(picturesBis)
    shaped_picturesBis = shaped_pictures
    
    # flattened pictures. flat_pictures is a list of numpy array with dimensions (nrPixels x 1)
    flat_pictures = flattenPictures(shaped_picturesBis)
    
    return labels, pictures, shaped_pictures, flat_pictures

