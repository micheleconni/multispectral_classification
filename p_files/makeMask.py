# -*- coding: utf-8 -*-
import nrrd
import numpy as np

# read file
filename = 'T02.nrrd'
maskname = 'mask.nrrd'

readdata, options = nrrd.read(filename)
imageShape = (readdata.shape)
x = imageShape[0]
y = imageShape[1]
z = imageShape[2]

# make mask
mask = np.zeros((x,y,z)).astype(int)
for i in range(0,x):
    for j in range(0,y):
        for k in range(0,z):
            mask[i][j][k]=1  

# write mask to nrrd file
nrrd.write(maskname, mask)

print ('Mask shape', mask.shape)
