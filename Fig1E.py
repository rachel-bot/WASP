# R. Brunetti, 210331
# Read in relative x positions of WASP and clathrin.
# Histogram the data and compare single cell means using a paired t-test.
# Figure 1 shows the panels used in paper Fig 1E
# p-value and its source single cell means are shown in Figure 2.

import os
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage import filters
import pandas as pd
import seaborn as sns
from scipy import stats
plt.ion()

nbins = 10 #number of bins in spatial histogram

#load cell contour
cell_contour = np.load('ex_cell_contour.npy')

#load data
data = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 1E')

#divide contour into regions
#set cell front to zero position
cell_contour[:,1] =  cell_contour[:,1]  - np.min(cell_contour[:,1]) * np.ones(np.shape(cell_contour[:,1]))
#set cell bottom to zero position
cell_contour[:,0] =  cell_contour[:,0]  - np.min(cell_contour[:,0]) * np.ones(np.shape(cell_contour[:,0]))
cell_front_x = np.min(cell_contour[:,1])
cell_rear_x = np.max(cell_contour[:,1])
spacing = cell_rear_x - cell_front_x
bins = [cell_front_x + i*spacing/nbins for i in range(nbins+1)]

#make a histogram of puncta x-coordinates
punctanums_means_wasp = np.histogram(data.WASPx, np.arange(0,1+1/nbins,1/(nbins)))[0]
total_puncta_wasp = np.sum(punctanums_means_wasp)
punctanums_means_CLTA = np.histogram(data.CLTAx, np.arange(0,1+1/nbins,1/(nbins)))[0]
total_puncta_CLTA = np.sum(punctanums_means_CLTA)

#make a binary filled cell
cellim = np.zeros((int(np.ceil(np.max(cell_contour[:,0]))),int(np.ceil(np.max(cell_contour[:,1])))))
for entry in cell_contour:
    cellim[int(entry[0]), int(entry[1])] = 1
cellim_binary = ndimage.binary_fill_holes(cellim).astype(int)
cellim_wasp = cellim_binary.astype(float)
cellim_CLTA = cellim_binary.astype(float)

#fill spatial bins with histogram values
for i in range(nbins-1):
    cellim_wasp[:, int(bins[i]):int(bins[i+1])] = cellim_wasp[:, int(bins[i]):int(bins[i+1])] * 100* punctanums_means_wasp[i]/total_puncta_wasp
    cellim_CLTA[:, int(bins[i]):int(bins[i+1])] = cellim_CLTA[:, int(bins[i]):int(bins[i+1])] * 100* punctanums_means_CLTA[i]/total_puncta_CLTA

#plot WASP histogram
plt.figure(figsize = (8,5))
ax1 = plt.subplot(2,2,1)
plt.title('WASP')
cmap_choice = plt.cm.get_cmap('YlGn')
cmap_choice.set_bad(color = 'white')
#smooth image
sigma = 7
blurred = ndimage.gaussian_filter1d(cellim_wasp, sigma)
blurred = blurred*cellim_binary + cellim_binary
blurred[blurred == 0] = np.NaN
blurred = blurred+cellim_binary
plt.imshow(blurred, cmap = cmap_choice, vmin = 0, vmax = 30)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
plt.axis('off')
ax1.invert_xaxis()

#plot CLTA histogram
ax2 = plt.subplot(2,2,3)
plt.title('CLTA')
cmap_choice = plt.cm.get_cmap('RdPu')
cmap_choice.set_bad(color = 'white')
blurred = ndimage.gaussian_filter1d(cellim_CLTA, sigma)
blurred = blurred*cellim_binary +cellim_binary
blurred[blurred == 0] = np.NaN
blurred = blurred+cellim_binary
plt.imshow(blurred, cmap = cmap_choice, vmin = 0, vmax = 30)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
plt.axis('off')
ax2.invert_xaxis()


#plot line scan
plt.subplot(2,2,(2,4))
xaxis = np.arange(0.05,1,1/(nbins))
plt.plot(xaxis, np.flip(100*punctanums_means_wasp/total_puncta_wasp), 'g-', linewidth = 2)
plt.plot(xaxis, np.flip(100*punctanums_means_CLTA/total_puncta_CLTA), 'm-', linewidth = 2)
plt.ylabel('percentage of puncta in bin')
plt.xlabel('relative distance from rear')
print('no. WASP puncta = ' + str(total_puncta_wasp))
print('no. CLTA puncta = ' + str(total_puncta_CLTA))


#plot all data points
plt.figure()
df_tot = pd.DataFrame({})
df_tot['x'] = np.concatenate([data.WASPx, data.CLTAx])
df_tot['ID'] = ['WASP' for i in data.WASPx] + ['CLTA' for i in data.CLTAx]
df_tot['x_flipped'] = np.ones(len(df_tot)) - df_tot.x
ax = sns.swarmplot(x = 'x_flipped', y = 'ID', data = df_tot, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1)

#isolate cell means
unique_cells = set(data.Cell)
cell_means_WASP = [1-np.nanmean(data[data.Cell == x].WASPx)for x in unique_cells]
cell_means_CLTA = [1-np.nanmean(data[data.Cell == x].CLTAx)for x in unique_cells]

#make a data frame for seaborn plotting
df_means = pd.DataFrame({})
df_means['x-position'] = np.concatenate([cell_means_WASP, cell_means_CLTA])
df_means['marker'] = ['WASP' for i in cell_means_WASP] + ['CLTA' for i in cell_means_CLTA]
df_means['cell_num'] = ['cell ' + str(i) for i in range(len(cell_means_WASP))] + ['cell ' + str(i) for i in range(len(cell_means_CLTA))]
g = sns.swarmplot(x = 'x-position', y = 'marker', hue = 'cell_num', data = df_means, size = 10, linewidth = 2, palette = 'cubehelix')
g.legend(bbox_to_anchor = (1.3, 0.5), loc = 7)
g.set_xlabel('Relative distance from cell rear')
plt.tight_layout()

#perform paired stats on cell means
#confirm there are enough observations per cell to do stats on cell level:
WASP_obs_per_cell = [np.sum(np.isfinite(data[data.Cell == x].WASPx))for x in unique_cells]
CLTA_obs_per_cell = [np.sum(np.isfinite(data[data.Cell == x].CLTAx))for x in unique_cells]
print('Average WASP puncta observed per cell = '+ str(int(np.mean(WASP_obs_per_cell))))
print('Average CLTA puncta observed per cell = '+ str(int(np.mean(CLTA_obs_per_cell))))

statistic, pvalue = stats.ttest_rel(cell_means_WASP, cell_means_CLTA)
print('pvalue = ' + str(pvalue))

x1, x2 = 0, 1
y, h, col = 1+ 0.05, 0.05, 'k'
plt.plot([y, y+h, y+h, y], [x1, x1, x2, x2], lw=1.5, c=col)
if pvalue > 0.001:
    plt.text(y+1.5*h, 0.685, "P = "+str("{:.3f}".format(pvalue)), ha='center', va='bottom', color=col, rotation = 270)
else:
    plt.text(y+1.5*h, 0.685, "P = " + str("{:.2e}".format(pvalue)), ha='center', va='bottom', color=col, rotation = 270)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#print cell means ± standard error
print('mean relative WASP position: ' + str(np.round(np.mean(cell_means_WASP), 2)) + ' ± ' + str(np.round(stats.sem(cell_means_WASP), 2)))
print('mean relative clathrin position: ' + str(np.round(np.mean(cell_means_CLTA), 2)) + ' ± ' + str(np.round(stats.sem(cell_means_CLTA), 2)))
