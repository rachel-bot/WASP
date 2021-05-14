# R. Brunetti, 210331
# Calculate the mean relative intensity of WASP signal at beads in each tenth
# of the cell. Underlying data is relative WASP signal at beads as a function
# of their distance from the leading edge, which was calculated for beads that
# the cell migrates the majority of its length over. Script produces paper
# figure 4B.
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure, morphology, filters
import pandas as pd
import seaborn as sns
import math
from scipy import stats, ndimage
import random
plt.ion()

# number of bins for histogram
nbins = 10
cell_contour = np.load('ex_cell_contour.npy')

# load in distance and relative intensity values
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 4B and S5AFig')

# put cell contour against x = 0 and y = 0
cell_contour[:,1] =  cell_contour[:,1]  - np.min(cell_contour[:,1]) * np.ones(np.shape(cell_contour[:,1]))
cell_contour[:,0] =  cell_contour[:,0]  - np.min(cell_contour[:,0]) * np.ones(np.shape(cell_contour[:,0]))

#divide contour into regions
cell_front_x = np.min(cell_contour[:,1])
cell_rear_x = np.max(cell_contour[:,1])
spacing = cell_rear_x - cell_front_x
cell_bins = [cell_front_x + i*spacing/nbins for i in range(nbins+1)]

#make the cell shaped binary to fill
cellim = np.zeros((int(np.ceil(np.max(cell_contour[:,0]))),int(np.ceil(np.max(cell_contour[:,1])))))
for entry in cell_contour:
    cellim[int(entry[0]), int(entry[1])] = 1
cellim_binary = ndimage.binary_fill_holes(cellim).astype(int)

#calculate mean relative WASP intensity in each 10% spatial bin
#function to discretize data
def round_nearest(x, a):
    return np.round(np.round(x / a) * a, -int(math.floor(math.log10(a))))

#discretize position to nearest tenth of the cell for the histogram
df['Position'] = round_nearest(df.relDist,0.1)

#calculate the intensity means in each region
bins = np.linspace(0,1,nbins+1)
mean_ints = []
for val in bins:
    val = np.round(val,2)
    mean_ints.append(df[df.Position == val].mean().relInt)


#fill cell regions with mean intensity values in each 10% of the cell
mean_ints = mean_ints/np.sum(mean_ints) #normalize for a density value
for i in range(nbins):
    cellim[:, int(cell_bins[i]):int(cell_bins[i+1])] = mean_ints[i]*100

#display spatial histogram
plt.figure(figsize = (8,5))
cmap_choice = plt.cm.get_cmap('magma')
cmap_choice.set_bad(color = 'white')
plt.imshow(cellim, cmap = cmap_choice)
cb = plt.colorbar(orientation = 'horizontal')
cb.set_ticks([])
plt.axis('off')


#smooth for nicer looking image and remove background
sigma = 7
blurred = ndimage.gaussian_filter1d(cellim, sigma)
blurred = blurred*cellim_binary +cellim_binary
blurred[blurred == 0] = np.NaN
blurred = blurred+cellim_binary
ax = plt.imshow(blurred, cmap = cmap_choice)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
plt.gca().invert_xaxis()
