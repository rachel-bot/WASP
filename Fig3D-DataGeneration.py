# R. Brunetti, 210331
# Script that allows user to select the left and right side of the invagination
# neck to track closing. Values are measured in pixels between selected edges.
# Values do not need to be exact, as we only want to identify approx. frame of closing.
# Otsu thresholding is also applied to determine relative intensity over time--
# particular interest in values at beginning of movie (open state, t_i) and at
# time of closing (t_f). Diameters and relative intensities are transferred to
# Excel for comparison  of intensities at t_i and t_f across multiple closing
# invaginations.
import os
from skimage import io, filters
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

#load time series of a single closing membrane invagination
imseq = io.imread('Fig3D_ex.tif')
membrane = imseq[:,0]
WASP = imseq[:,1]

#for each frame, click the left and right of the bud.
#if closed, double click in the same position, making diameter (x1-x0) = 0
i = 0
bud_diameter = []
for image in membrane:
    plt.imshow(image, cmap = 'Greys_r')
    i += 1
    plt.title('Click the left and right of the neck; frame = ' + str(i) +'/'+str(np.shape(membrane)[0]))
    bud_edges = plt.ginput(2)
    x0 = bud_edges[0][0]
    x1 = bud_edges[1][0]
    bud_diameter.append(x1-x0)
    plt.close()

print('neck diameter: ')
for x in bud_diameter:
    print(np.round(x,2))


# Otsu threshold WASP based on intensities from all frames
thresh = filters.threshold_multiotsu(np.concatenate(WASP))
binaryWASP = []
for image in WASP:
    image = ndimage.gaussian_filter(image,1)
    image[image<thresh[-1]] = 0
    image[image!=0] = 1
    binaryWASP.append(image)

#mask WASP and calculate total intensity of non-background signal
masked_WASP = binaryWASP*WASP
int_intensity_masked = [np.sum(x) for x in masked_WASP]
#normalize for each cell to compare across cells
int_intensity_masked_norm = int_intensity_masked/np.max(int_intensity_masked)

print('WASP intensity values: ')
for x in int_intensity_masked_norm:
    print(np.round(x,2))
