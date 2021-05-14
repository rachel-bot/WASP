# R. Brunetti, 210331
# Script that compares WASP puncta appearance in WT and Cdc42 KO cells
# Due to limited observations per cell, statistics were applied to replicate means.
# Position of appearance was not significantly different in these cell backgrounds
# via an unpaired two-tailed t-test. Figure here is paper Fig 4F.
import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import ndimage, stats
import matplotlib.pyplot as plt
plt.ion()

rep1color = '#E79D06'
rep2color = '#53B7E8'
rep3color = '#009F72'

# Load data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 4F')

# separate cell lines
df_WT = df[df.Condition == 'WT']
df_KO = df[df.Condition == 'KO']

# separate data out by replicate
dataWT1 = df_WT[df_WT.Replicate == 1]
dataWT2 = df_WT[df_WT.Replicate == 2]
dataWT3 = df_WT[df_WT.Replicate == 3]
dataWT = pd.concat([dataWT1, dataWT2, dataWT3])

dataKO1 = df_KO[df_KO.Replicate == 1]
dataKO2 = df_KO[df_KO.Replicate == 2]
dataKO3 = df_KO[df_KO.Replicate == 3]
dataKO = pd.concat([dataKO1, dataKO2, dataKO3])

# adding fields for seaborn plotting
dataWT['ID'] = ['WT' for i in range(len(dataWT))]
dataWT['Replicate'] = ['1' for i in range(len(dataWT1))]+['2' for i in range(len(dataWT2))]+['3' for i in range(len(dataWT3))]
dataKO['ID'] = ['Cdc42 KO' for i in range(len(dataKO))]
dataKO['Replicate'] = ['1' for i in range(len(dataKO1))]+['2' for i in range(len(dataKO2))]+['3' for i in range(len(dataKO3))]
datatot = pd.concat([dataWT, dataKO])

#swapping x axis to be distance from cell rear
f = plt.figure()
datatot['inverse_x'] = np.ones(np.shape(datatot.x)) - datatot.x
ax = sns.swarmplot(x = 'inverse_x', y = 'ID', data = datatot, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1)
plt.xlabel('distance from cell rear')
plt.ylabel('Condition')

#isolate cell means. Are there enough observations per cell to support stats on cell means?
#no, there were not.
unique_cells_WT_rep1 = set(dataWT1.cell)
WT_obs_per_cell_1 = [len(dataWT1[dataWT1.cell == x]) for x in unique_cells_WT_rep1]
WT_singlecellmean_1 = [np.mean(dataWT1[dataWT1.cell == x].x) for x in unique_cells_WT_rep1]
unique_cells_WT_rep2 = set(dataWT2.cell)
WT_obs_per_cell_2 = [len(dataWT2[dataWT2.cell == x]) for x in unique_cells_WT_rep2]
WT_singlecellmean_2 = [np.mean(dataWT2[dataWT2.cell == x].x) for x in unique_cells_WT_rep2]
unique_cells_WT_rep3 = set(dataWT3.cell)
WT_obs_per_cell_3 = [len(dataWT3[dataWT3.cell == x]) for x in unique_cells_WT_rep3]
WT_singlecellmean_3 = [np.mean(dataWT3[dataWT3.cell == x].x) for x in unique_cells_WT_rep3]
WT_singlecellmeans = np.concatenate([WT_singlecellmean_1, WT_singlecellmean_2, WT_singlecellmean_3])

unique_cells_KO_rep1 = set(dataKO1.cell)
KO_obs_per_cell_1 = [len(dataKO1[dataKO1.cell == x]) for x in unique_cells_KO_rep1]
KO_singlecellmean_1 = [np.mean(dataKO1[dataKO1.cell == x].x) for x in unique_cells_KO_rep1]
unique_cells_KO_rep2 = set(dataKO2.cell)
KO_obs_per_cell_2 = [len(dataKO2[dataKO2.cell == x]) for x in unique_cells_KO_rep2]
KO_singlecellmean_2 = [np.mean(dataKO2[dataKO2.cell == x].x) for x in unique_cells_KO_rep2]
unique_cells_KO_rep3 = set(dataKO3.cell)
KO_obs_per_cell_3 = [len(dataKO3[dataKO3.cell == x]) for x in unique_cells_KO_rep3]
KO_singlecellmean_3 = [np.mean(dataKO3[dataKO3.cell == x].x) for x in unique_cells_KO_rep3]
KO_singlecellmeans = np.concatenate([KO_singlecellmean_1, KO_singlecellmean_2, KO_singlecellmean_3])

print('Average puncta observed per WT cell = '+ str(int(np.mean(np.concatenate([WT_obs_per_cell_1, WT_obs_per_cell_2, WT_obs_per_cell_3])))))
print('Average puncta observed per KO cell = '+ str(int(np.mean(np.concatenate([KO_obs_per_cell_1, KO_obs_per_cell_2, KO_obs_per_cell_3])))))
print('Not enough observations per cell to do stats on single cells. Doing on replicate means instead.')

#calculate replicate means (account for flipping xaxis in plotting)
WTmeans = [1-dataWT1.x.mean(), 1-dataWT2.x.mean(), 1-dataWT3.x.mean()]
KOmeans = [1-dataKO1.x.mean(), 1-dataKO2.x.mean(), 1-dataKO3.x.mean()]

#perform unpaired t-test on replicate means
statistic, pvalue = stats.ttest_ind(WTmeans, KOmeans)
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

#plot means
plt.scatter(WTmeans[0], -0.07, color = rep1color, linewidth = 1, edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(WTmeans[1], 0, color = rep2color,  linewidth = 1, edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(WTmeans[2], 0.07, color = rep3color,  linewidth = 1, edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(KOmeans[0], 0.93, color = rep1color,  linewidth = 1, edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(KOmeans[1], 1, color = rep2color,  linewidth = 1, edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(KOmeans[2], 1.07, color = rep3color,  linewidth = 1, edgecolor = 'k', zorder = 200, s = 75)

plt.xlabel('Relative distance from cell rear')
plt.suptitle('Position of puncta appearance')
plt.tight_layout()

#print means and standard errors
print('Mean position of appearance in WT cells: ' + str(np.round(np.mean(WTmeans),2)) + ' ± ' + str(np.round(stats.sem(WTmeans),2)))
print('Mean position of appearance in KO cells: ' + str(np.round(np.mean(KOmeans),2)) + ' ± ' + str(np.round(stats.sem(KOmeans),2)))
