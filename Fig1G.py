# R. Brunetti, 210331
# Script that separates out +/- EDTA puncta position disappearance
# and compares single cells means via an unpaired t-test
# Left side of Figure 1 is used in paper Fig1G
# and statistics shown in Figure 2 are also used in Fig1G

import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import ndimage, stats
import matplotlib.pyplot as plt
plt.ion()

Replicate1color = '#E79D06'
Replicate2color = '#53B7E8'
Replicate3color = '#009F72'
mypal = {'replicate 1':Replicate1color, 'replicate 2':Replicate2color, 'replicate 3':Replicate3color}

nbins = 10 #dividing cell into 10% bins but can be changed

# Load data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 1G')

# separate +/- EDTA conditions
df_EDTA = df[df.Condition == 'EDTA']
df_noEDTA = df[df.Condition == 'no EDTA']


# separate data out by condition
dataEDTA1 = df_EDTA[df_EDTA.Replicate == 1]
dataEDTA2 = df_EDTA[df_EDTA.Replicate == 2]
dataEDTA3 = df_EDTA[df_EDTA.Replicate == 3]
dataEDTA = pd.concat([dataEDTA1, dataEDTA2, dataEDTA3])


datanoEDTA1 = df_noEDTA[df_noEDTA.Replicate == 1]
datanoEDTA2 = df_noEDTA[df_noEDTA.Replicate == 2]
datanoEDTA3 = df_noEDTA[df_noEDTA.Replicate == 3]
datanoEDTA = pd.concat([datanoEDTA1, datanoEDTA2, datanoEDTA3])

# adding fields for seaborn plotting
dataEDTA['ID'] = ['+EDTA' for i in range(len(dataEDTA))]
dataEDTA['Replicate'] = ['1' for i in range(len(dataEDTA1))]+['2' for i in range(len(dataEDTA2))]+['3' for i in range(len(dataEDTA3))]
datanoEDTA['ID'] = ['-EDTA' for i in range(len(datanoEDTA))]
datanoEDTA['Replicate'] = ['1' for i in range(len(datanoEDTA1))]+['2' for i in range(len(datanoEDTA2))]+['3' for i in range(len(datanoEDTA3))]
datatot = pd.concat([dataEDTA, datanoEDTA])

# load in cell contour for cell-shaped histogram
cell_contour = np.load('ex_cell_contour.npy')
#set cell front to zero position
cell_contour[:,1] =  cell_contour[:,1]  - np.min(cell_contour[:,1]) * np.ones(np.shape(cell_contour[:,1]))
#set cell bottom to zero position
cell_contour[:,0] =  cell_contour[:,0]  - np.min(cell_contour[:,0]) * np.ones(np.shape(cell_contour[:,0]))

# divide puncta positions into bins
bins0 = np.linspace(0,1,nbins + 1)

gEDTA = np.histogram(dataEDTA.x, bins = bins0)
nEDTA = gEDTA[0]/np.sum(gEDTA[0])

gnoEDTA = np.histogram(datanoEDTA.x, bins = bins0)
nnoEDTA = gnoEDTA[0]/np.sum(gnoEDTA[0])


#divide contour into regions
cell_front_x = np.min(cell_contour[:,1])
cell_rear_x = np.max(cell_contour[:,1])

spacing = cell_rear_x - cell_front_x
bins = [cell_front_x + i*spacing/nbins for i in range(nbins+1)]

#readjusting coords from anywhere in space to a bounded matrix
cellim = np.zeros((int(np.ceil(np.max(cell_contour[:,0]))),int(np.ceil(np.max(cell_contour[:,1])))))
for entry in cell_contour:
    cellim[int(entry[0]), int(entry[1])] = 1

cellim_binary = ndimage.binary_fill_holes(cellim).astype(int)
cellim = cellim_binary.astype(float)

cellimEDTA = np.copy(cellim)
cellimnoEDTA = np.copy(cellim)

# fill in cell-shaped histogram with densities for each marker
for i in range(nbins-1):
    cellimEDTA[:, int(bins[i]):int(bins[i+1])] = cellimEDTA[:, int(bins[i]):int(bins[i+1])] * 100* gEDTA[0][i]/np.sum(gEDTA[0])
    cellimnoEDTA[:, int(bins[i]):int(bins[i+1])] = cellimnoEDTA[:, int(bins[i]):int(bins[i+1])] * 100* gnoEDTA[0][i]/np.sum(gnoEDTA[0])

# display histograms
plt.figure(figsize = (8,5))
ax1 = plt.subplot(2,2,1)
cmap_choice = plt.cm.get_cmap('magma')
cmap_choice.set_bad(color = 'white')
plt.imshow(cellimEDTA, cmap = cmap_choice)
plt.axis('off')
plt.title('+EDTA')
sigma = 7
blurred = ndimage.gaussian_filter1d(cellimEDTA, sigma)
blurred = blurred*cellim_binary +cellim_binary
blurred[blurred == 0] = np.NaN
blurred = blurred+cellim_binary
plt.imshow(blurred, cmap = cmap_choice)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 1)
plt.imshow(blurred, cmap = cmap_choice, vmin = 0, vmax = 25)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
plt.axis('off')
ax1.invert_xaxis()


ax2 = plt.subplot(2,2,3)
cmap_choice = plt.cm.get_cmap('magma')
cmap_choice.set_bad(color = 'white')
plt.imshow(cellimnoEDTA, cmap = cmap_choice)
plt.axis('off')
plt.title('-EDTA')
sigma = 7
blurred = ndimage.gaussian_filter1d(cellimnoEDTA, sigma)
blurred = blurred*cellim_binary +cellim_binary
blurred[blurred == 0] = np.NaN
blurred = blurred+cellim_binary
plt.imshow(blurred, cmap = cmap_choice)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 1)
plt.imshow(blurred, cmap = cmap_choice, vmin = 0, vmax = 25)
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
plt.axis('off')
ax2.invert_xaxis()


# line scan
ax3 = plt.subplot(2,2,(2,4))
xaxis = np.linspace(0,1,nbins+1)
xaxis = xaxis+0.05
plt.plot(np.flip(xaxis[:-1]), 100*gEDTA[0]/np.sum(gEDTA[0]), '-', linewidth = 2, label = '+EDTA')
plt.plot(np.flip(xaxis[:-1]), 100*gnoEDTA[0]/np.sum(gnoEDTA[0]), '-', linewidth = 2, label = '-EDTA')
plt.legend()
plt.ylabel('Percentage of puncta in bin')
plt.xlabel('Relative distance from cell rear')


#swapping x axis to be distance from cell rear
f = plt.figure()
datatot['inverse_x'] = np.ones(np.shape(datatot.x)) - datatot.x
ax = sns.swarmplot(x = 'inverse_x', y = 'ID', data = datatot, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1)
plt.xlabel('distance from cell rear')
plt.ylabel('Condition')

#isolate cell means. Are there enough observations per cell to support stats on cell means? yes
unique_cells_EDTA_rep1 = set(dataEDTA1.cell)
EDTA_obs_per_cell_1 = [len(dataEDTA1[dataEDTA1.cell == x]) for x in unique_cells_EDTA_rep1]
EDTA_singlecellmean_1 = [np.mean(dataEDTA1[dataEDTA1.cell == x].x) for x in unique_cells_EDTA_rep1]
unique_cells_EDTA_rep2 = set(dataEDTA2.cell)
EDTA_obs_per_cell_2 = [len(dataEDTA2[dataEDTA2.cell == x]) for x in unique_cells_EDTA_rep2]
EDTA_singlecellmean_2 = [np.mean(dataEDTA2[dataEDTA2.cell == x].x) for x in unique_cells_EDTA_rep2]
unique_cells_EDTA_rep3 = set(dataEDTA3.cell)
EDTA_obs_per_cell_3 = [len(dataEDTA3[dataEDTA3.cell == x]) for x in unique_cells_EDTA_rep3]
EDTA_singlecellmean_3 = [np.mean(dataEDTA3[dataEDTA3.cell == x].x) for x in unique_cells_EDTA_rep3]
EDTA_singlecellmeans = np.concatenate([EDTA_singlecellmean_1, EDTA_singlecellmean_2, EDTA_singlecellmean_3])

unique_cells_noEDTA_rep1 = set(datanoEDTA1.cell)
noEDTA_obs_per_cell_1 = [len(datanoEDTA1[datanoEDTA1.cell == x]) for x in unique_cells_noEDTA_rep1]
noEDTA_singlecellmean_1 = [np.mean(datanoEDTA1[datanoEDTA1.cell == x].x) for x in unique_cells_noEDTA_rep1]
unique_cells_noEDTA_rep2 = set(datanoEDTA2.cell)
noEDTA_obs_per_cell_2 = [len(datanoEDTA2[datanoEDTA2.cell == x]) for x in unique_cells_noEDTA_rep2]
noEDTA_singlecellmean_2 = [np.mean(datanoEDTA2[datanoEDTA2.cell == x].x) for x in unique_cells_noEDTA_rep2]
unique_cells_noEDTA_rep3 = set(datanoEDTA3.cell)
noEDTA_obs_per_cell_3 = [len(datanoEDTA3[datanoEDTA3.cell == x]) for x in unique_cells_noEDTA_rep3]
noEDTA_singlecellmean_3 = [np.mean(datanoEDTA3[datanoEDTA3.cell == x].x) for x in unique_cells_noEDTA_rep3]
noEDTA_singlecellmeans = np.concatenate([noEDTA_singlecellmean_1, noEDTA_singlecellmean_2, noEDTA_singlecellmean_3])

print('Average puncta observed per +EDTA cell = '+ str(int(np.mean([EDTA_obs_per_cell_1, EDTA_obs_per_cell_2, EDTA_obs_per_cell_3]))))
print('Average puncta observed per -EDTA cell = '+ str(int(np.mean([noEDTA_obs_per_cell_1, noEDTA_obs_per_cell_2, noEDTA_obs_per_cell_3]))))

#perform unpaired t-test on single cell means
statistic, pvalue = stats.ttest_ind(EDTA_singlecellmeans, noEDTA_singlecellmeans)
print('pvalue = ' + str(pvalue))

x1, x2 = 0, 1
y, h, col = 1+ 0.05, 0.05, 'k'
plt.plot([y, y+h, y+h, y], [x1, x1, x2, x2], lw=1.5, c=col)
if pvalue > 0.001:
    plt.text(y+1.5*h, 0.685, "P = "+str(float(round(pvalue, 3))), ha='center', va='bottom', color=col, rotation = 270)
else:
    plt.text(y+1.5*h, 0.685, "P = " + str("{:.2e}".format(pvalue)), ha='center', va='bottom', color=col, rotation = 270)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#make a data frame for seaborn plotting of means
df_means = pd.DataFrame({})
df_means['x-position'] = np.concatenate([1-EDTA_singlecellmeans, 1-noEDTA_singlecellmeans])
df_means['condition'] = ['EDTA' for i in EDTA_singlecellmeans] + ['noEDTA' for i in noEDTA_singlecellmeans]
df_means['replicate'] = ['replicate 1' for i in EDTA_singlecellmean_1] + ['replicate 2' for i in EDTA_singlecellmean_1] + ['replicate 3' for i in EDTA_singlecellmean_1] + ['replicate 1' for i in noEDTA_singlecellmean_1] + ['replicate 2' for i in noEDTA_singlecellmean_1] + ['replicate 3' for i in noEDTA_singlecellmean_1]
g = sns.swarmplot(x = 'x-position', y = 'condition', hue = 'replicate', data = df_means, size = 10, linewidth = 2, palette = mypal)
g.legend(bbox_to_anchor = (1.3, 0.5), loc = 7, borderaxespad = -3)
g.set_xlabel('Relative distance from cell rear')
plt.suptitle('Position of puncta disappearance')
plt.tight_layout()

#print means and standard error
print('realtive position of disappearance w/o EDTA: ' + str(np.round(np.mean(1-noEDTA_singlecellmeans),2)) + ' ± ' + str(np.round(stats.sem(noEDTA_singlecellmeans),2)))
print('realtive position of disappearance w/ EDTA: ' + str(np.round(np.mean(1-EDTA_singlecellmeans),2)) + ' ± ' + str(np.round(stats.sem(EDTA_singlecellmeans),2)))
