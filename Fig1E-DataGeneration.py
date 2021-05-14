# R. Brunetti, 210331
# Identify WASP and clathrin puncta from a difference of Gaussians image.
# Extract puncta positions and calculate their relative position along the cell length.
# Print data to combine across cells and analyze.

#Use example image, which is equal to sample Cell 1, Replicate 1 in
#Excel sheet Fig 1E

import os
import skimage
from skimage import io, filters, measure, morphology
import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
from scipy import ndimage
plt.ion()

####USER SPECIFIED VALUES####
cellmask_thresh_scale = 0.3 #scale cell mask segmentation to make sure it isolates cell correctly
num_sections = 10 #will be doing region classification based on intervals of 10% but can be changed

#load data
fname = 'Fig1E_ex.tif'

#import and split single cell image into channels
ss_im = skimage.io.imread(fname)
CLTA = ss_im[0]
WASP = ss_im[1]


#use blurred WASP signal to make a cell mask
WASP_max = ndimage.maximum_filter(WASP, size = 5, mode = 'constant')
WASP_blurred = filters.gaussian(WASP_max, sigma = 5)
plt.subplot(1,2,1)
plt.imshow(WASP_blurred)
thresh = filters.threshold_multiotsu(WASP_blurred)
binary = WASP_blurred > cellmask_thresh_scale*thresh[1]
plt.subplot(1,2,2)
plt.imshow(binary)
plt.suptitle('Assess cell mask. Press enter to continue')
input()
plt.close()

#keep only the largest region
membrane_binary = np.copy(binary)
labels = skimage.measure.label(membrane_binary)
regions = skimage.measure.regionprops(labels)
cell_index = np.argmax([x.area for x in regions])  + 1
cell_mask = np.zeros(np.shape(membrane_binary))
cell_mask[labels == cell_index] = 1
label_cell_mask = skimage.measure.label(cell_mask)
rp = skimage.measure.regionprops(label_cell_mask)
val = rp[np.argmax([x.area for x in rp])].label
cell_mask[label_cell_mask != val] = 0
plt.imshow(cell_mask)

#make contour of cell_mask
contours = skimage.measure.find_contours(cell_mask, 0.8)
COI = np.argmax(([len(x) for x in contours]))
cell_contour = contours[COI]
plt.plot(cell_contour[:,1], cell_contour[:,0])

#find the centroid
cell_centroid = regions[cell_index-1].centroid
plt.plot(cell_centroid[1], cell_centroid[0], '.')

#find the front and the back x-coordinates of the cell
#x is second column
cell_front_x = np.min(cell_contour[:,1])
cell_rear_x = np.max(cell_contour[:,1])
front_index = np.argmin(cell_contour[:,1])
rear_index = np.argmax(cell_contour[:,1])

plt.plot(cell_contour[front_index][1], cell_contour[front_index][0], 'r.')
plt.plot(cell_contour[rear_index][1], cell_contour[rear_index][0], 'm.')
plt.suptitle('Assess bounds and centroid. Press enter to continue.')
input()
plt.close()

#split the cell into sections
spacing_x = cell_rear_x - cell_front_x
bins_x = [cell_front_x + i*spacing_x/num_sections for i in range(num_sections+1)]

#binary of WASP puncta
WASP_binary = skimage.filters.difference_of_gaussians(np.copy(WASP), 1, 2)
WASP_binary = WASP_binary*cell_mask
WASPthresh = skimage.filters.threshold_multiotsu(WASP_binary)
WASP_binary[WASP_binary<WASPthresh[1]] = 0
WASP_binary[np.nonzero(WASP_binary)] = 1

#repeat for CLTA
CLTA_binary = skimage.filters.difference_of_gaussians(np.copy(CLTA), 1, 2)
CLTA_binary = CLTA_binary*cell_mask
CLTAthresh = skimage.filters.threshold_multiotsu(CLTA_binary)
CLTA_binary[CLTA_binary<CLTAthresh[1]] = 0
CLTA_binary[np.nonzero(CLTA_binary)] = 1

#determine centroids of puncta
labels_WASP = skimage.measure.label(WASP_binary)
WASP_puncta_props = skimage.measure.regionprops(labels_WASP)
WASP_puncta_centroids = [x.centroid for x in WASP_puncta_props]
labels_CLTA = skimage.measure.label(CLTA_binary)
CLTA_puncta_props = skimage.measure.regionprops(labels_CLTA)
CLTA_puncta_centroids = [x.centroid for x in CLTA_puncta_props]

#Let's see how we did:
plt.figure()
plt.subplot(2,3,1)
plt.title('WASP')
plt.imshow(WASP)
plt.axis('off')

plt.subplot(2,3,2)
plt.title('WASP + puncta labels')
plt.imshow(WASP_binary)
for entry in WASP_puncta_centroids:
    plt.plot(entry[1], entry[0], 'ro', lw = 2, fillstyle = 'none')
plt.axis('off')


#map xs onto the cell axis
cell_vector = [cell_contour[rear_index][1]-cell_contour[front_index][1], cell_contour[rear_index][0]-cell_contour[front_index][0]]
x_vals = []
for centroid in WASP_puncta_centroids:
    puncta_vector = [centroid[1]-cell_contour[front_index][1], centroid[0]-cell_contour[front_index][0]]
    fractional_x = np.dot(cell_vector, puncta_vector)/np.linalg.norm(cell_vector)/np.linalg.norm(cell_vector)
    x_vals.append(fractional_x)

ax1 = plt.subplot(2,3,3)
plt.title('WASP distribution')
plt.hist(x_vals)
plt.xlim([0,1])
ax1.set_aspect(0.03)
plt.tight_layout()

#to transfer to a csv
print('WASP coords')
for x in x_vals:
    print(x)

plt.subplot(2,3,4)
plt.title('CLTA')
plt.imshow(CLTA)
plt.axis('off')

plt.subplot(2,3,5)
plt.title('CLTA + puncta labels')
plt.imshow(CLTA_binary)
for entry in CLTA_puncta_centroids:
    plt.plot(entry[1], entry[0], 'ro', lw = 2, fillstyle = 'none')
plt.axis('off')

#map xs onto the cell axis
x_vals_CLTA = []
for centroid in CLTA_puncta_centroids:
    puncta_vector = [centroid[1]-cell_contour[front_index][1], centroid[0]-cell_contour[front_index][0]]
    fractional_x = np.dot(cell_vector, puncta_vector)/np.linalg.norm(cell_vector)/np.linalg.norm(cell_vector)
    x_vals_CLTA.append(fractional_x)

ax2 = plt.subplot(2,3,6)
plt.title('CLTA distributions')
plt.hist(x_vals_CLTA)
plt.xlim([0,1])
ax2.set_aspect(0.03)
plt.tight_layout()

#to transfer to a csv
print('CLTA coords')
for x in x_vals_CLTA:
    print(x)

#put histograms on same axis
max_ax1 = ax1.get_ylim()[1]
max_ax2 = ax2.get_ylim()[1]
if max_ax1 > max_ax2:
    ax2.set_ylim([0, max_ax1])
else:
    ax1.set_ylim([1, max_ax2])
