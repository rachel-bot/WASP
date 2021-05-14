# R. Brunetti, 210331
# Script to measure relative WASP signal as a function of distance from the
# leading edge. User selects a bead that traverses the cell length (determined
# in FIJI). Frames where this bead is under the cell are isolated. For each
# of these frames, user specifies the cell front and cell rear which are used
# to calculate the cell major axis. The bead coordinates are projected onto this
# line to give distance from the leading edge. At the same time, the WASP signal
# in the vicinity of the selected bead is isolated and the integrated intensity
# is calculated for each time point and normalized to the max frame. Data is
# printed for transfer to Excel.

# Example image corresponds to cell 4, bead 1 in Excel sheet "Fig 4B and S5A Fig"
# as there is only one bead that passes under the cell during this movie. Note
# that segmentation here is not perfect, but all we care about is extracting
# the cell front and cell rear.
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure, morphology, filters
import math
plt.ion()

#user defined values
boxradius = 3 #area to isolate around the bead
cellmask_thresh_scale = 0.5 #scaling factor for making cell mask
FOI = 5 #frame of interest you would like to select bead in
#if watch segmentation is true, you can determine how well thethresholding
#with the specified scaling is performing
watchSegmentation = True


#load imaging data
fname = 'Fig4B_ex'
imseq = io.imread(fname + '.tif')
WASP = imseq[:,1]
beads = imseq[:,0]

#display the frame of interest and select a single bead that traverses the whole
#the cell (determined through watching movies in FIJI)
FOI = FOI - 1
RGB = np.zeros((np.shape(imseq)[2], np.shape(imseq)[3],3))
RGB[:,:,0] = WASP[FOI]/np.max(WASP[FOI])
RGB[:,:,1] = 5*beads[FOI]/np.max(5*beads[FOI])

#select bead of interest for this analysis
plt.figure(figsize = (12,7))
plt.imshow(RGB)
plt.title('select the bead of interest')
plt.axis('off')
plt.tight_layout()
bead_centroid = plt.ginput(1)
plt.close()

#find all the frames where the bead is under the cell
cellmask = []
for frame in WASP:
    wasp_blurred = filters.gaussian(np.copy(frame), sigma = 5)
    thresh = filters.threshold_multiotsu(wasp_blurred)
    binary = wasp_blurred > cellmask_thresh_scale*thresh[1]
    if watchSegmentation:
            plt.subplot(1,2,1)
            plt.imshow(wasp_blurred)
            plt.axis('off')
            plt.title('blurred WASP signal')
            plt.subplot(122)
            plt.imshow(binary)
            plt.axis('off')
            plt.title('cell mask')
            plt.pause(0.2)
            plt.close()
    cellmask.append(binary)

#For each frame the bead is under the cell, select the cell rear and then cell front
#for relative length calculations
i = -1
coords = []
frame_present = []
for frame in cellmask:
    i+=1
    if frame[int(bead_centroid[0][1]), int(bead_centroid[0][0])] == 1:
        plt.imshow(frame)
        plt.axis('off')
        plt.title('Click the cell rear and then the cell front')
        plt.plot(int(bead_centroid[0][0]), int(bead_centroid[0][1]), 'r.')
        plt.draw()
        coords.append(plt.ginput(2))
        frame_present.append(i)

#project coordinate onto cell major axis
distance_along_axis = []
for coord_set in coords:
    AB = coord_set
    #vector from rear to front
    v1 = [AB[0][0]- AB[1][0], AB[0][1]-AB[1][1]]
    #vector from bead to front
    v2 = [bead_centroid[0][0]- AB[1][0], bead_centroid[0][1]-AB[1][1]]
    distance_along_axis.append(np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v1))


#isolate WASP and bead signal at selection over time
bead_ROI = beads[frame_present,int(bead_centroid[0][1]-boxradius):int(bead_centroid[0][1]+boxradius), int(bead_centroid[0][0]-boxradius):int(bead_centroid[0][0]+boxradius)]
WASP_ROI = WASP[frame_present,int(bead_centroid[0][1]-boxradius):int(bead_centroid[0][1]+boxradius), int(bead_centroid[0][0]-boxradius):int(bead_centroid[0][0]+boxradius)]
WASP_int = [np.sum(x) for x in WASP_ROI]

#display normalized WASP signal as a fucntion of distance from the cell front
plt.figure()
plt.subplot(1,2,1)
plt.plot(distance_along_axis, WASP_int/np.max(WASP_int))
plt.xlabel('Distance from the cell front')
plt.ylabel('Normalized WASP signal')
plt.subplot(1,2,2)
plt.imshow(np.concatenate(WASP_ROI))
plt.axis('off')
plt.xlabel('WASP signal at bead')

#print values to transfer to Excel
#return relative position and WASP intensity
print('relative distance from cell front')
for i in range(len(distance_along_axis)):
    print(distance_along_axis[i])

print('relative WASP intensity')
for i in range(len(WASP_int)):
    print(WASP_int[i]/np.max(WASP_int))
