# -*- coding: utf-8 -*-

""" A script to make an excel-file with the paths to the tiff files """

import numpy as np
#import nrrd
#import scipy.io as sio
import os 
import os.path
import pandas as pd
import csv


""" Settings """
numberOfImages = 5


""" Make dataframe and excel file """

length = np.zeros(shape = (numberOfImages,3))
tabell = pd.DataFrame(data = length, columns = ['ID', 'Image', 'Mask']) # Create the dataframe
directory = os.getcwd() # Set the current directory
j = 0
for i in range(0,numberOfImages+1):
    maskName = '\mask.nrrd' #Making a name for the mask
    if i < 10:
        fileName = 'T0' + str(i) + '.nrrd' # Making a file name to check if the file is in the folder
        imageName = '\T0' + str(i) + '.nrrd'# Making a name for the new image file
    else:
        fileName = 'T' + str(i) + '.nrrd'
        imageName = '\T' + str(i) + '.nrrd'
    if os.path.isfile(fileName):
        tabell.iloc[j,0] = i        # A column for the image ID
        tabell.iloc[j,1] = directory +  imageName # A colummn for the path to the images
        tabell.iloc[j,2] = directory +  maskName  # A column for the path to the mask of the image
        j += 1  
    print(i)
    
tabell.to_csv('paths.csv') # Make a excel-file with the paths. 
