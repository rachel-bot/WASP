# R. Brunetti, 210331
# Calculate the Pearson correlation coefficient for
# a bead and the protein of interest in a 1.5x1.5 um
# ROI around the bead. Record output values in a csv for
# comparison across different proteins of interest.

#Example image corresponds to the first WASP entry in Excel sheet Fig 2G
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

#load zoomed image of protein of interest at a bead
os.chdir('Fig2G')
imseq = skimage.io.imread('Fig2G_CLTA_ex.tif')
bead = imseq[1]
protein_of_interest = imseq[0]

# flatten and normalize both channels
nonzero_protein_of_interest = [x for x in protein_of_interest.ravel() if x != 0]
nonzero_protein_of_interest = nonzero_protein_of_interest/np.max(nonzero_protein_of_interest)
nonzero_bead = [x for x in bead.ravel() if x != 0]
nonzero_bead = nonzero_bead/np.max(nonzero_bead)

# rename variables for standardized variable names
x = nonzero_bead
y = nonzero_protein_of_interest

#fit linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
def linefitline(m, b, x):
    return b*np.ones(np.shape(x)) + m*np.ones(np.shape(x)) * x
line = linefitline(slope, intercept, nonzero_bead)

#plot the pixel kernel density estimate and fitted line
xy = np.vstack([x,y])
z = scipy.stats.gaussian_kde(xy)(xy)
plt.scatter(x, y, c=z, s=100, edgecolor='')
plt.plot(nonzero_bead,line,'r:')

#return the correlation coefficient to record
print('Pearson correlation coefficient: ' + str(r_value))
