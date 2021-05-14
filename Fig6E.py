# R. Brunetti, 210331
# Display and compare the ratio of distance traveled along nanoridged patterns
# (x) and perpendicular to nanoridged patterns (y) across cell backgrounds
#(WT v WASP KO lines). Makes paper Fig 6E.

import pickle as pkl
import os
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import seaborn as sns
import numpy as np
plt.ion()

#scale to convert pixels to um (imaged using a 20x objective)
umperpx = 0.65

#load data and scale to um
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 6D and E')
df['Distance'] = df['Distance']*umperpx
df['x Distance'] = df['x Distance']*umperpx
df['y Distance'] = df['y Distance']*umperpx
# add ratio field
df['ratio_xy'] = df['x Distance']/df['y Distance']


#separate out replicates and cell lines to easily calculate means
rep1 = df[df.replicate == 1]
WT_rep1 = rep1[rep1['Cell Type'] == 'Wild Type']
KO_rep1 = rep1[rep1['Cell Type'] == 'WASP KO']
KO2_rep1 = rep1[rep1['Cell Type'] == 'WASP KO 2']

rep2 = df[df.replicate == 2]
WT_rep2 = rep2[rep2['Cell Type'] == 'Wild Type']
KO_rep2 = rep2[rep2['Cell Type'] == 'WASP KO']
KO2_rep2 = rep2[rep2['Cell Type'] == 'WASP KO 2']

rep3 = df[df.replicate == 3]
WT_rep3 = rep3[rep3['Cell Type'] == 'Wild Type']
KO_rep3 = rep3[rep3['Cell Type'] == 'WASP KO']
KO2_rep3 = rep3[rep3['Cell Type'] == 'WASP KO 2']

#isolate cell line means
WT_xy_means = [WT_rep1.ratio_xy.mean(), WT_rep2.ratio_xy.mean(), WT_rep3.ratio_xy.mean()]
KO_xy_means = [KO_rep1.ratio_xy.mean(), KO_rep2.ratio_xy.mean(), KO_rep3.ratio_xy.mean()]
KO2_xy_means = [KO2_rep1.ratio_xy.mean(), KO2_rep2.ratio_xy.mean(), KO2_rep3.ratio_xy.mean()]


#plot
sns.swarmplot(x = 'Cell Type', y = 'ratio_xy', data = df, size=10, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1, alpha = 1)

plt.scatter(-0.15, WT_xy_means[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 200)
plt.scatter(0, WT_xy_means[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 200)
plt.scatter(0.15, WT_xy_means[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 200)

plt.scatter(0.85, KO_xy_means[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 200)
plt.scatter(1, KO_xy_means[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 200)
plt.scatter(1.15, KO_xy_means[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 200)

plt.scatter(1.85, KO2_xy_means[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 200)
plt.scatter(2, KO2_xy_means[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 200)
plt.scatter(2.15, KO2_xy_means[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 200)

x1, x2 = 0, 1
y, h, col = 16, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

x1, x2 = 0, 2
y, h, col = 17.5, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

###statistics###
#paired t-test on replicate means
#WT x/y ratio v KO1/2 x/y ratio
s,p_xy_WT_KO = stats.ttest_rel(WT_xy_means, KO_xy_means)
s,p_xy_WT_KO2 = stats.ttest_rel(WT_xy_means, KO2_xy_means)
print('p xy ratio between WT v KO: ' + str(p_xy_WT_KO))
print('p xy ratio between WT v KO2: ' + str(p_xy_WT_KO2))

#print means and standard error
print('WT mean: ')
print(str(np.round(np.mean(WT_xy_means),2)) + ' ± ' + str(np.round(stats.sem(WT_xy_means),2)))
print('KO mean: ')
print(str(np.round(np.mean(KO_xy_means),2)) + ' ± ' + str(np.round(stats.sem(KO_xy_means),2)))
print('KO2 mean: ')
print(str(np.round(np.mean(KO2_xy_means),2)) + ' ± ' + str(np.round(stats.sem(KO2_xy_means),2)))
