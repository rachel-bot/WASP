# R. Brunetti, 210331
# Separate the cell length into 10% bins. Calculate frequency of the
# x-coordinates of WASP and Arp3 puncta along the cell length for each bin in
# co-expressing cells. Display distribution as a cell-shaped histogram. Compare
# the mean position of WASP and Arp3 puncta at the cell level using a paired
# two-tailed t-test. Generates paper figure 4C and the p-value reported.
import os
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import ndimage
from skimage import filters
import pandas as pd
import seaborn as sns
from scipy import stats
plt.ion()

#number of bins to divide cell into
nbins = 10

#load in and reposition the cell outline
cell_contour = np.load('ex_cell_contour.npy')
#set cell front to zero position
cell_contour[:,1] =  cell_contour[:,1]  - np.min(cell_contour[:,1]) * np.ones(np.shape(cell_contour[:,1]))
#set cell bottom to zero position
cell_contour[:,0] =  cell_contour[:,0]  - np.min(cell_contour[:,0]) * np.ones(np.shape(cell_contour[:,0]))

#load WASP and Arp3 x-coordinate data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 5C')

#divide contour into regions
cell_front_x = np.min(cell_contour[:,1])
cell_rear_x = np.max(cell_contour[:,1])
spacing = cell_rear_x - cell_front_x
bins = [cell_front_x + i*spacing/nbins for i in range(nbins+1)]

#calculate number of puncta per bin and normalize by total counts to make a density
punctanums_means_wasp = np.histogram(df['x-wasp'], np.arange(0,1,1/(nbins+1)))[0]
total_puncta_wasp = np.sum(punctanums_means_wasp)
punctanums_means_arp = np.histogram(df['x-arp'], np.arange(0,1,1/(nbins+1)))[0]
total_puncta_arp = np.sum(punctanums_means_arp)

#make the binary cell shape
cellim = np.zeros((int(np.ceil(np.max(cell_contour[:,0]))),int(np.ceil(np.max(cell_contour[:,1])))))
for entry in cell_contour:
    cellim[int(entry[0]), int(entry[1])] = 1
cellim_binary = ndimage.binary_fill_holes(cellim).astype(int)

#map histogram values onto cell shape for each channel
cellim_wasp = cellim_binary.astype(float)
cellim_arp = cellim_binary.astype(float)
for i in range(nbins-1):
    cellim_wasp[:, int(bins[i]):int(bins[i+1])] = cellim_wasp[:, int(bins[i]):int(bins[i+1])] * 100* punctanums_means_wasp[i]/total_puncta_wasp
    cellim_arp[:, int(bins[i]):int(bins[i+1])] = cellim_arp[:, int(bins[i]):int(bins[i+1])] * 100* punctanums_means_arp[i]/total_puncta_arp

#display histograms
plt.figure(figsize = (8,5))
plt.subplot(2,2,1)
cmap_choice = plt.cm.get_cmap('magma')
cmap_choice.set_bad(color = 'white')
plt.imshow(cellim_wasp, cmap = cmap_choice)
plt.title('WASP')
cb = plt.colorbar(orientation = 'horizontal')
cb.set_ticks([])
plt.axis('off')
sigma = 7
blurred = scipy.ndimage.gaussian_filter1d(cellim_wasp, sigma)
blurred = blurred*cellim_binary +cellim_binary
blurred[blurred == 0] = np.NaN
blurred = blurred+cellim_binary
plt.imshow(blurred, cmap = cmap_choice, vmin = 0, vmax = 25)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
plt.axis('off')
plt.gca().invert_xaxis()

plt.subplot(2,2,3)
cmap_choice = plt.cm.get_cmap('magma')
cmap_choice.set_bad(color = 'white')
plt.imshow(cellim_arp, cmap = cmap_choice)
plt.title('Arp3')
cb = plt.colorbar(orientation = 'horizontal')
cb.set_ticks([])
plt.axis('off')
blurred = scipy.ndimage.gaussian_filter1d(cellim_arp, sigma)
blurred = blurred*cellim_binary +cellim_binary
blurred[blurred == 0] = np.NaN
blurred = blurred+cellim_binary
plt.imshow(blurred, cmap = cmap_choice, vmin = 0, vmax = 25)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
plt.axis('off')
plt.gca().invert_xaxis()

#display a cdf
plt.subplot(2,2,2)
cdf_wasp = np.cumsum(punctanums_means_wasp/total_puncta_wasp)
cdf_arp = np.cumsum(punctanums_means_arp/total_puncta_arp)
xaxis = np.arange(0,1,1/(nbins))
plt.plot(xaxis, np.flip(100*punctanums_means_wasp/total_puncta_wasp), 'g-', linewidth = 2)
plt.plot(xaxis, np.flip(100*punctanums_means_arp/total_puncta_arp), 'm-', linewidth = 2)
plt.ylabel('percentage of puncta in bin')

#display a line scan across the histogram
plt.subplot(2,2,4)
plt.plot(xaxis, cdf_wasp, 'g-', linewidth = 2)
plt.plot(xaxis, cdf_arp, 'm-', linewidth = 2)
plt.ylabel('probability')
plt.xlabel('relative distance from rear')
plt.subplots_adjust(wspace = 0.4)
print('no. WASP puncta = ' + str(total_puncta_wasp))
print('no. Arp3 puncta = ' + str(total_puncta_arp))

#second figure showing all x data points
plt.figure()
df_tot = pd.DataFrame({})
df_tot['x'] = np.concatenate([df['x-wasp'], df['x-arp']])
df_tot['ID'] = ['WASP' for i in range(len(df['x-wasp']))] + ['Arp3' for i in range(len(df['x-arp']))]
df_tot['x_flipped'] = np.ones(len(df_tot)) - df_tot.x
ax = sns.swarmplot(x = 'x_flipped', y = 'ID', data = df_tot, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1)
statistic, pvalue = stats.mannwhitneyu(df['x-wasp'], df['x-arp'])

#Are there enough observations per cell to warrant stats on single cells?
#yep!
df_1 = df[df.replicate == 1]
df_2 = df[df.replicate == 2]
df_3 = df[df.replicate == 3]
#WASP observations per cell
wasp_obs_per_cell_rep1 = [np.sum(np.isfinite(df_1[df_1.field == x]['x-wasp'])) for x in np.unique(df_1.field)]
wasp_obs_per_cell_rep2 = [np.sum(np.isfinite(df_2[df_2.field == x]['x-wasp'])) for x in np.unique(df_2.field)]
wasp_obs_per_cell_rep3 = [np.sum(np.isfinite(df_3[df_3.field == x]['x-wasp'])) for x in np.unique(df_3.field)]
print('WASP observations per cell:')
print(wasp_obs_per_cell_rep1 + wasp_obs_per_cell_rep2 + wasp_obs_per_cell_rep3)
#Arp3 observations per cell
arp_obs_per_cell_rep1 = [np.sum(np.isfinite(df_1[df_1.field == x]['x-arp'])) for x in np.unique(df_1.field)]
arp_obs_per_cell_rep2 = [np.sum(np.isfinite(df_2[df_2.field == x]['x-arp'])) for x in np.unique(df_2.field)]
arp_obs_per_cell_rep3 = [np.sum(np.isfinite(df_3[df_3.field == x]['x-arp'])) for x in np.unique(df_3.field)]
print('arp observations per cell:')
print(arp_obs_per_cell_rep1 + arp_obs_per_cell_rep2 + arp_obs_per_cell_rep3)

#mean value per marker per cell
mean_wasp_x_per_cell_rep1 = [np.nanmean(df_1[df_1.field == x]['x-wasp']) for x in np.unique(df_1.field)]
mean_wasp_x_per_cell_rep2 = [np.nanmean(df_2[df_2.field == x]['x-wasp']) for x in np.unique(df_2.field)]
mean_wasp_x_per_cell_rep3 = [np.nanmean(df_3[df_3.field == x]['x-wasp']) for x in np.unique(df_3.field)]
mean_wasp_x = mean_wasp_x_per_cell_rep1 + mean_wasp_x_per_cell_rep2 + mean_wasp_x_per_cell_rep3

mean_arp_x_per_cell_rep1 = [np.nanmean(df_1[df_1.field == x]['x-arp']) for x in np.unique(df_1.field)]
mean_arp_x_per_cell_rep2 = [np.nanmean(df_2[df_2.field == x]['x-arp']) for x in np.unique(df_2.field)]
mean_arp_x_per_cell_rep3 = [np.nanmean(df_3[df_3.field == x]['x-arp']) for x in np.unique(df_3.field)]
mean_arp_x = mean_arp_x_per_cell_rep1 + mean_arp_x_per_cell_rep2 + mean_arp_x_per_cell_rep3

#plot means
df_means = pd.DataFrame({})
df_means['x-position'] = np.concatenate([1-np.asarray(mean_wasp_x), 1-np.asarray(mean_arp_x)])
df_means['marker'] = ['WASP' for i in mean_wasp_x] + ['Arp3' for i in mean_arp_x]
df_means['cell_num'] = ['cell ' + str(i) for i in range(len(mean_wasp_x))] + ['cell ' + str(i) for i in range(len(mean_arp_x))]
g = sns.swarmplot(x = 'x-position', y = 'marker', hue = 'cell_num', data = df_means, size = 10, linewidth = 2, palette = 'cubehelix')
g.legend(bbox_to_anchor = (1.3, 0.5), loc = 7)
g.set_xlabel('Relative distance from cell rear')
plt.tight_layout()

#paired t-test on single cell means
s, pvalue = stats.ttest_rel(mean_wasp_x, mean_arp_x)
print('pvalue = ' + str(pvalue))

#display stats value
x1, x2 = 0, 1
y, h, col = 1+ 0.05, 0.05, 'k'
plt.plot([y, y+h, y+h, y], [x1, x1, x2, x2], lw=1.5, c=col)
plt.text(y+1.5*h, 0.685, "P = "+str(float(round(pvalue, 3))), ha='center', va='bottom', color=col, rotation = 270)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#display means and standard error
print('Mean relative WASP position: ' + str(1-np.round(np.mean(mean_wasp_x),2)) + ' ± ' + str(np.round(stats.sem(mean_wasp_x),2)))
print('Mean relative Arp3 position: ' + str(1-np.round(np.mean(mean_arp_x),2)) + ' ± ' + str(np.round(stats.sem(mean_arp_x),2)))
