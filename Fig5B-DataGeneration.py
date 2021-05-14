#R. Brunetti, 210331
#Calculate Pearson correlation coefficients between WASP and Arp3.
#Generate a negative control by calculating the coefficient between
#WASP and a 90 degree rotation of Arp3.

#Example data corresponds to replicate 1, ROI 1 which is line 2 for the
#regular correlation and line 22 for the rotated correlation 
import os
import skimage
from skimage import io, filters
import numpy as np
import math
import scipy
import csv
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats
import random
plt.ion()

# load two marker ROI
imseq = skimage.io.imread('Fig5B_ex.tif')
wasp = imseq[0]
arp = imseq[1]

# flatten and normalize both channels
nonzero_arp = [x for x in arp.ravel() if x != 0]
nonzero_arp = nonzero_arp/np.max(nonzero_arp)
nonzero_WASP = [x for x in wasp.ravel() if x != 0]
nonzero_WASP = nonzero_WASP/np.max(nonzero_WASP)

# rename variables with standardized variable names
x = nonzero_WASP
y = nonzero_arp

#fit linear regression
slope, intercept, r_value1, p_value, std_err = stats.linregress(x,y)
def linefitline(m, b, x):
    return b*np.ones(np.shape(x)) + m*np.ones(np.shape(x)) * x
line1 = linefitline(slope, intercept, nonzero_WASP)

#plot the pixel kernel density estimate and fitted line
plt.subplot(1,2,1)
xy = np.vstack([x,y])
z = scipy.stats.gaussian_kde(xy)(xy)
plt.scatter(x, y, c=z, s=100, edgecolor='')
plt.plot(nonzero_WASP,line1,'r:')
plt.title('WASP-Arp correlation')

#create a negative control where we rotate one channel 90 degrees
m = np.copy(arp)
#rotating 90 degrees
arp_reshuffle = [[m[j][i] for j in range(len(m))] for i in range(len(m[0])-1,-1,-1)]
#redefining y to be the transformed image (x is the same)
y =  [x for x in np.asarray(arp_reshuffle).ravel() if x != 0]
y = y/np.max(y)

#fit the linear regression between original x and transformed y
slope, intercept, r_value2, p_value, std_err = stats.linregress(nonzero_WASP, y)
line2 = linefitline(slope, intercept, nonzero_WASP)

#plot the pixel kernel density estimate and fitted line for the control
plt.subplot(1,2,2)
xy = np.vstack([x,y])
z = scipy.stats.gaussian_kde(xy)(xy)
plt.scatter(x, y, c=z, s=100, edgecolor='')
plt.plot(nonzero_WASP,line2, 'r:')
plt.title('negative control; Arp rotated')

print('Pearson correlation coefficient: ' + str(r_value1))
print('Pearson correlation coefficient of negative control: ' + str(r_value2))
