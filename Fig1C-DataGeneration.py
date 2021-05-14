# R. Brunetti, 210331
# Script that categorizes WASP puncta as the leading edge (LE), cell front,
# cell rear, or uropod (U) based on the following definitions:
# LE = area eroded by a 20X20 structuring element in the front 40% of the cell mask
# cell front = front 50% of the cell front excluding the LE
# cell rear = 50% to 80% of the cell length
# uropod = rear 20% of the cell
# classifications were chosen for robustness across cells

# variables with "region" = WASP signal in that region
# variable with "footprint" = cell mask (for area measure)

#Use example image, which is equal to sample image_prefix = '0823-002' in
#Excel sheet Fig 1C

import os
import skimage
from skimage import io, filters, measure, morphology
import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import pickle as pkl
plt.ion()


####USER SPECIFIED VALUES####
cell_outline_thresh = 0.28 #scale cell mask segmentation to make sure it isolates cell correctly
num_sections = 10 #will be doing region classification based on intervals of 10%, but can be changed

#load microscopy data
#cell should be oriented left to right
fname = 'Fig1C_ex.tif'

#import and split single cell image into channels
ss_im = skimage.io.imread(fname)
membrane = ss_im[1]
WASP = ss_im[1]

#threshold membrane for binarized cell mask
thresh_mem = 0.3*np.max(membrane) #starting guess for threshold
membrane_binary = np.copy(membrane)
membrane_binary[membrane_binary < cell_outline_thresh*thresh_mem] = 0 #scaling threshold by user input "cell_outline_thresh"
membrane_binary[np.nonzero(membrane_binary)] = 1
plt.subplot(1,2,1)
plt.imshow(membrane)
plt.subplot(1,2,2)
plt.imshow(membrane_binary)
plt.suptitle('Assess thresholding and press enter to continue')
input()
plt.close()


#keep largest object in the binary and make a mask
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

#find the centroid of the cell
cell_centroid = regions[cell_index-1].centroid
plt.plot(cell_centroid[1], cell_centroid[0], '.')

#find the front and the back x-coordiantes of the cell
#x is second column
cell_front_x = np.min(cell_contour[:,1])
cell_rear_x = np.max(cell_contour[:,1])
front_index = np.argmin(cell_contour[:,1])
rear_index = np.argmax(cell_contour[:,1])

plt.plot(cell_contour[front_index][1], cell_contour[front_index][0], 'r.')
plt.plot(cell_contour[rear_index][1], cell_contour[rear_index][0], 'm.')
plt.suptitle('Assess cell front, back and centroid. Press enter to continue')
input()
plt.close()

#split the cell into sections along the major axis
spacing_x = cell_rear_x - cell_front_x
bins_x = [cell_front_x + i*spacing_x/num_sections for i in range(num_sections+1)]

#make a WASP DoG and threshold using Otsu thresholding
WASP_copy = np.copy(WASP)
WASP_DoG = filters.difference_of_gaussians(WASP_copy, 1, 2)
thresh = filters.threshold_multiotsu(WASP_DoG)
WASP_DoG[WASP_DoG > thresh[-1]] = 1
WASP_DoG[WASP_DoG != 1] = 0
WASP = WASP * WASP_DoG
masked_WASP = cell_mask*WASP #WASP puncta contained in the cell mask

#calculate leading edge (LE) enrichment by eroding the front 40% of the cell
forty = bins_x[4]
kernel = 20 #determined empirically
mask_eroded =  ndimage.morphology.binary_erosion(cell_mask, structure = np.ones((kernel, kernel)))
periphery = masked_WASP - mask_eroded*WASP #WASP puncta in cell outline
periphery_footprint = cell_mask*WASP_copy - mask_eroded*WASP_copy #raw WASP signal at cell outline
LE_footprint = periphery_footprint[:, 0:int(forty)] #raw WASP signal at cell outline in front 40%
LE_region = periphery[:, 0:int(forty)] #WASP puncta at cell outline in front 40%

#remove LE WASP puncta from the rest of the signal
periphery_binary = np.copy(LE_region)
periphery_binary[periphery_binary>0] = 1
periphery_binary = np.ones(np.shape(periphery_binary))-periphery_binary
masked_WASP_no_LE = masked_WASP[:,0:int(forty)]*periphery_binary
masked_WASP_no_LE = np.hstack([masked_WASP_no_LE, masked_WASP[:,int(forty):]])

#remove LE region from the rest of the footprint
periphery_binary_footprint = np.copy(LE_footprint)
periphery_binary_footprint[periphery_binary_footprint>0] = 1
periphery_binary_footprint = np.ones(np.shape(periphery_binary_footprint))-periphery_binary_footprint
masked_WASP_footprint = cell_mask*WASP_copy
masked_WASP_no_LE_footprint = masked_WASP_footprint[:,0:int(forty)]*periphery_binary_footprint
masked_WASP_no_LE_footprint = np.hstack([masked_WASP_no_LE_footprint, masked_WASP_footprint[:,int(forty):]])

#identify cell front (front 50% excluding LE)
fifty = bins_x[5]
front_region = masked_WASP_no_LE[:,:int(fifty)]
front_region_footprint = masked_WASP_no_LE_footprint[:,:int(fifty)]

#make figure and plot cell front and use it to set intensity scale for other regions
fig = plt.figure()
plt.subplot(2,4,2)
plt.imshow(front_region_footprint, cmap = 'gray')
vmin, vmax = plt.gci().get_clim()
vmax = 0.4*vmax
plt.imshow(front_region_footprint, vmin = vmin, vmax = vmax, cmap = 'gray')
plt.axis('off')
plt.title('Cell Front')
ax2 = plt.subplot(2,4,6)
plt.imshow(front_region, vmin = vmin, vmax = 0.1*vmax, cmap = 'gray')
plt.axis('off')
front_signal = np.sum(front_region)
front_norm = front_signal/np.count_nonzero(front_region_footprint)

#plot LE
ax1 = plt.subplot(2,4,5)
plt.imshow(LE_region, vmin = vmin, vmax = 0.1*vmax, cmap = 'gray')
plt.axis('off')
plt.subplot(2,4,1)
plt.imshow(LE_footprint, vmin = vmin, vmax = vmax, cmap = 'gray')
plt.axis('off')
plt.title('Leading Edge')
LE_signal = np.sum(LE_region)
LE_norm = LE_signal/np.count_nonzero(LE_footprint)

#identify cell rear (between 50 and 80%) and plot
eighty = bins_x[8]
rear_region = masked_WASP_no_LE[:,int(fifty):int(eighty)]
rear_region_footprint = masked_WASP_no_LE_footprint[:,int(fifty):int(eighty)]
ax3 = plt.subplot(2,4,7)
plt.imshow(rear_region, vmin = vmin, vmax = 0.1*vmax, cmap = 'gray')
plt.axis('off')
plt.subplot(2,4,3)
plt.imshow(rear_region_footprint, vmin = vmin, vmax = vmax, cmap = 'gray')
plt.axis('off')
plt.title('Cell Rear')
rear_signal = np.sum(rear_region)
rear_norm = rear_signal/np.count_nonzero(rear_region_footprint)

#identify uropod (rear 20%) and plot
U_region = masked_WASP_no_LE[:,int(eighty):]
U_region_footprint = masked_WASP_no_LE_footprint[:,int(eighty):]
ax4 = plt.subplot(2,4,8)
plt.imshow(U_region, vmin = vmin, vmax = 0.1*vmax, cmap = 'gray')
plt.axis('off')
plt.subplot(2,4,4)
plt.imshow(U_region_footprint, vmin = vmin, vmax = vmax, cmap = 'gray')
plt.axis('off')
plt.title('Uropod')
U_signal = np.sum(U_region)
U_norm = U_signal/np.count_nonzero(U_region_footprint)

#make regions binary to get number of pixels ("area")
front_bin = front_region>0
LE_bin = LE_region>0
rear_bin = rear_region>0
U_bin = U_region>0

#calculate region properties of labeled binary to get puncta centroids for plotting
front_rp = measure.regionprops(measure.label(front_bin*1))
for entry in front_rp:
    centroid = entry.centroid
    ax2.plot(centroid[1], centroid[0], 'ro', lw = 2, fillstyle = 'none')
LE_rp = measure.regionprops(measure.label(LE_bin*1))
for entry in LE_rp:
    centroid = entry.centroid
    ax1.plot(centroid[1], centroid[0], 'ro', lw = 2, fillstyle = 'none')
rear_rp = measure.regionprops(measure.label(rear_bin*1))
for entry in rear_rp:
    centroid = entry.centroid
    ax3.plot(centroid[1], centroid[0], 'ro', lw = 2, fillstyle = 'none')
U_rp = measure.regionprops(measure.label(U_bin*1))
for entry in U_rp:
    centroid = entry.centroid
    ax4.plot(centroid[1], centroid[0], 'ro', lw = 2, fillstyle = 'none')

#calculate percent signal contribution of each region
tot_signal = LE_signal + front_signal + rear_signal + U_signal
print('fraction LE:' + str(np.round(LE_signal/tot_signal,2)))
print('fraction front:' + str(np.round(front_signal/tot_signal,2)))
print('fraction rear:' + str(np.round(rear_signal/tot_signal,2)))
print('fraction U:' + str(np.round(U_signal/tot_signal,2)))

#calculate area normalized signal (percent signal / percent size of cell)
tot_size = np.count_nonzero(LE_footprint) + np.count_nonzero(front_region_footprint) + np.count_nonzero(rear_region_footprint) + np.count_nonzero(U_region_footprint)
LE_norm = LE_signal/tot_signal * tot_size/np.count_nonzero(LE_footprint)
front_norm = front_signal/tot_signal * tot_size/np.count_nonzero(front_region_footprint)
rear_norm = rear_signal/tot_signal * tot_size/np.count_nonzero(rear_region_footprint)
U_norm = U_signal/tot_signal * tot_size/np.count_nonzero(U_region_footprint)

#calculate what percentage each of these is relative to the total summed area-normalized signal
new_base = LE_norm+front_norm+rear_norm+U_norm
LE_norm2 = np.round(LE_norm/new_base, 2)
front_norm2 = np.round(front_norm/new_base, 2)
rear_norm2 = np.round(rear_norm/new_base, 2)
U_norm2 = np.round(U_norm/new_base, 2)

print('norm LE: ' + str(LE_norm2))
print('norm front: ' + str(front_norm2))
print('norm rear: ' + str(rear_norm2))
print('norm U: ' + str(U_norm2))
