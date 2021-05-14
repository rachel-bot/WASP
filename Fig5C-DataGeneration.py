# R. Brunetti, 210331
# Script that segments the cell based on blurred WASP signal. User-adjusted
# Otsu thresholding is used to best identify puncta from the background in a
# difference of Gaussians image. x coordinates are extractdd. All cells are
# oriented left to right and the x coordinates are mapped between 0 (leftmost
# point) and 1 (rightmost point) along the cell legnth

#Example image corresponds to replicate 2, field 2 in Excel sheet Fig 5C

import os
import skimage
from skimage import io, filters, measure, morphology
import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import pickle as pkl
plt.ion()

#user defined values
#scales used to confirm good segmentation
cellmask_thresh_scale = 0.75
arpthresh_scale = 3.5
waspthresh_scale = 0.75
#number of sections to divide the cell into
num_sections = 10

#load data
fname = 'Fig5C_ex.tif'

#import and split single cell image
arp = skimage.io.imread(fname)[0]
wasp = skimage.io.imread(fname)[1]


# use wasp signal to make a cell mask
wasp_blurred = filters.gaussian(wasp, sigma = 5)
plt.subplot(1,2,1)
plt.imshow(wasp_blurred)
thresh = filters.threshold_multiotsu(wasp_blurred)
binary = wasp_blurred > cellmask_thresh_scale*thresh[1]
plt.subplot(122)
plt.imshow(binary)
plt.suptitle('Press enter is segmentation is okay. Else rescale.')
input()
plt.close()

#continue refining the cell mask
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

#find the front and the back of the cell
#x is second column
cell_front_x = np.min(cell_contour[:,1])
cell_rear_x = np.max(cell_contour[:,1])
front_index = np.argmin(cell_contour[:,1])
rear_index = np.argmax(cell_contour[:,1])

#display mask and properties and determine if accurate
plt.plot(cell_contour[front_index][1], cell_contour[front_index][0], 'r.')
plt.plot(cell_contour[rear_index][1], cell_contour[rear_index][0], 'm.')
plt.suptitle('Confirm that centroid and x coordinates of front and back look okay')
input()
plt.close()

#split the cell into sections
spacing_x = cell_rear_x - cell_front_x
bins_x = [cell_front_x + i*spacing_x/num_sections for i in range(num_sections+1)]

#Identify arp puncta
arp_binary = skimage.filters.difference_of_gaussians(np.copy(arp), 1, 2)
arp_binary = arp_binary*cell_mask
arpthresh = skimage.filters.threshold_multiotsu(arp_binary)
arp_binary[arp_binary<arpthresh_scale*arpthresh[1]] = 0
arp_binary[np.nonzero(arp_binary)] = 1

#determine puncta centroids
labels_arp = skimage.measure.label(arp_binary)
arp_puncta_props = skimage.measure.regionprops(labels_arp)
puncta_centroids = [x.centroid for x in arp_puncta_props]

#plot Arp data
plt.figure(figsize = (12,7))
plt.subplot(2,4,1)
plt.title('Arp3')
plt.axis('off')
plt.imshow(arp)

plt.subplot(2,4,2)
plt.title('Arp3 puncta')
plt.imshow(arp*cell_mask)
plt.axis('off')
for entry in puncta_centroids:
    plt.plot(entry[1], entry[0], 'r+')

plt.subplot(2,4,3)
plt.title('Thresholded Arp3 signal')
plt.imshow(arp_binary)
plt.axis('off')

#calculate and append relative distance of Arp3 puncta along cell length
cell_vector = [cell_contour[rear_index][1]-cell_contour[front_index][1], cell_contour[rear_index][0]-cell_contour[front_index][0]]
x_vals = []
for centroid in puncta_centroids:
    puncta_vector = [centroid[1]-cell_contour[front_index][1], centroid[0]-cell_contour[front_index][0]]
    fractional_x = np.dot(cell_vector, puncta_vector)/np.linalg.norm(cell_vector)/np.linalg.norm(cell_vector)
    x_vals.append(fractional_x)

#display distribution as histogram
ax1 = plt.subplot(2,4,4)
bins = np.linspace(0, 1, num_sections+1)
plt.hist(x_vals, bins = bins)

#print values to transfer to Excel
print('Arp3 puncta x-coords:')
for x in x_vals:
    print(x)

#repeat process for WASP. Indentify WASP puncta
wasp_binary = skimage.filters.difference_of_gaussians(np.copy(wasp), 1, 2)
wasp_binary = wasp_binary*cell_mask
waspthresh = skimage.filters.threshold_multiotsu(wasp_binary)
wasp_binary[wasp_binary<waspthresh_scale*waspthresh[1]] = 0
wasp_binary[np.nonzero(wasp_binary)] = 1

#determine puncta centroids
labels_wasp = skimage.measure.label(wasp_binary)
wasp_puncta_props = skimage.measure.regionprops(labels_wasp)
puncta_centroids = [x.centroid for x in wasp_puncta_props]

#plot WASP data
plt.subplot(2,4,5)
plt.imshow(wasp)
plt.axis('off')
plt.title('WASP')

plt.subplot(2,4,6)
plt.imshow(wasp*cell_mask)
plt.title('WASP puncta')
plt.axis('off')
for entry in puncta_centroids:
    plt.plot(entry[1], entry[0], 'r+')

plt.subplot(2,4,7)
plt.title('Thresholded WASP signal')
plt.imshow(wasp_binary)
plt.axis('off')

#calculate and append relative distance of WASP puncta along cell length
cell_vector = [cell_contour[rear_index][1]-cell_contour[front_index][1], cell_contour[rear_index][0]-cell_contour[front_index][0]]
x_vals = []
for centroid in puncta_centroids:
    puncta_vector = [centroid[1]-cell_contour[front_index][1], centroid[0]-cell_contour[front_index][0]]
    fractional_x = np.dot(cell_vector, puncta_vector)/np.linalg.norm(cell_vector)/np.linalg.norm(cell_vector)
    x_vals.append(fractional_x)

#display distribution as histogram
ax2 = plt.subplot(2,4,8)
plt.hist(x_vals, bins = bins)

#put histograms on same scale
if ax1.get_ylim()[1] > ax2.get_ylim()[1]:
    ax2.set_ylim([0, ax1.get_ylim()[1]])
else:
    ax1.set_ylim([0, ax2.get_ylim()[1]])

#print values to transfer to Excel
print('WASP puncta x-coords:')
for x in x_vals:
    print(x)
