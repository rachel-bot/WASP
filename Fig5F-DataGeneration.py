# R. Brunetti, 210331
# Script to measure relative background subtracted integrated intensity of Arp3
# signal at beads that move under the cell lamellipod. User specifies beads
# that will pass under a cell lamellipod. The region is tracked over time and
# and the frame where the highest value pixel is, is isolated. For this position
# the bead signal is used as a mask to separate signal at the bead from the
# background. The mean of the background signal is then subtracted off of each
# signal pixel and these values are summed for a background subtracted
# integrated intensity at the bead. Negative values sugest local depletion at
# the bead. We find that negative values and those less than ~three thousand AU
# correspond with noise as opposed to true Arp3 puncta regcruitment (as reflected
# in the percentages in Fig 5E).

#Example WT and KO movies included in Fig5EandF folder

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure, morphology, filters
import scipy
from scipy import optimize
import math
plt.ion()

#user defined values
boxradius = 5 #area to isolate around the bead
FOI = 15 #frame of interest you would like to select beads in
sizestep = 200 #dimension of zoomed inzet for easier bead selection

#load image
os.chdir('Fig5EandF')
fname = 'Fig5EandF_WT_ex_cropped.tif'

#order = frame, channel, x, y
imseq = io.imread(fname)
beads = imseq[:,1]
Arp3 = imseq[:,0]


#display the frame of interest
FOI = FOI - 1
RGB = np.zeros((np.shape(imseq)[2], np.shape(imseq)[3],3))
RGB[:,:,0] = Arp3[FOI]/np.max(Arp3[FOI])
RGB[:,:,1] = 5*beads[FOI]/np.max(5*beads[FOI])

#select the top left corner of the region to zoom in on
plt.figure(figsize = (12,7))
plt.imshow(RGB)
plt.title('choose left hand corner of zoom')
plt.axis('off')
plt.tight_layout()
plt.draw()

#specify number of beads you would like to analyze in this run
#keep number <= 6 so each puncta and signal kymograph can be evaluated
corner = plt.ginput(1)
plt.imshow(RGB[int(corner[0][1]):int(corner[0][1]+sizestep), int(corner[0][0]):int(corner[0][0]+sizestep),:])
plt.title('How many beads?')
plt.draw()
numbeads = input('Enter number of beads')

#select the beads
plt.draw()
bead_centroid = plt.ginput(int(numbeads))
plt.close()

#identify frame with max Arp3 pixel value for each bead
Arp3_ROI_tot = []
bead_ROI_tot = []
Arp3_int = []
Arp3_ROI_all_time = []
for i in range(len(bead_centroid)):
    Arp3_ROI = Arp3[:,int(bead_centroid[i][1]+corner[0][1]-boxradius):int(bead_centroid[i][1]+corner[0][1]+boxradius), int(bead_centroid[i][0]+corner[0][0]-boxradius):int(bead_centroid[i][0]+corner[0][0]+boxradius)]
    bead_ROI = beads[:,int(bead_centroid[i][1]+corner[0][1]-boxradius):int(bead_centroid[i][1]+corner[0][1]+boxradius), int(bead_centroid[i][0]+corner[0][0]-boxradius):int(bead_centroid[i][0]+corner[0][0]+boxradius)]
    Arp3_int.append(np.max([np.max(x) for x in Arp3_ROI]))
    index = np.argmax([np.max(x) for x in Arp3_ROI])
    Arp3_ROI_tot.append(Arp3_ROI[index])
    bead_ROI_tot.append(bead_ROI[index])
    Arp3_ROI_all_time.append(np.concatenate(Arp3_ROI))

#Otsu threshold bead signal
bead_ROI_tot_bin = []
for bead in bead_ROI_tot:
    thresh = filters.threshold_otsu(bead)
    bead_bin = np.copy(bead)
    bead_bin[bead_bin < thresh] = 0
    bead_bin[bead_bin != 0] = 1
    bead_ROI_tot_bin.append(bead_bin)


#calculate background and do background subtraction
intensity = []
Arp3_ROI_tot_bs = []
for i in range(len(Arp3_ROI_tot)):
    image = Arp3_ROI_tot[i]
    background = image - bead_ROI_tot_bin[i]*image
    average_background = np.mean(background[background != 0])
    foreground = bead_ROI_tot_bin[i]*image
    background_sub_int = np.sum(foreground) - average_background*np.sum(foreground != 0)
    intensity.append(background_sub_int)
    Arp3_ROI_tot_bs.append(foreground - average_background*np.ones(np.shape(foreground)))

#display Arp3 kymographs for each bead. Allows user to assess whether isolating
#the frame with the max Arp3 pixel did a good job of identifying the frame where
#the Arp3 puncta is brightest/most pronounced above background.
#We also tried identifying the frame based on integrated intensity but found
#this approach to be more robust.
plt.figure()
plt.title('Arp3 kymograph of signal at selected beads')
plt.imshow(np.concatenate(Arp3_ROI_all_time, axis = 1))
plt.axis('off')

#display Arp3 at the point where the integrated intensity is brightest and the
#background subtracted Arp3 signal at bead is to be reported
maxbeads = 6
plt.figure()
for i in range(len(Arp3_ROI_tot)):
    plt.subplot(2,maxbeads,i+1)
    plt.imshow(Arp3_ROI_tot[i])
    plt.axis('off')
    plt.subplot(2,maxbeads,i+1+maxbeads)
    plt.imshow(Arp3_ROI_tot_bs[i])
    plt.axis('off')
plt.suptitle('top = raw signal and bottom = background sub')


#print values for transfer to Excel
print('Background subtracted Arp3 integrated intensity:')
for value in intensity:
    print(int(value))
