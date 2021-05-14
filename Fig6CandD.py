# R. Brunetti, 210331
# Display and compare distance traveled in 2D and x and y across different
# cell backgrounds.
# Reading in Excel sheet "Fig 6C" makes paper figure 6C while reading in Excel
# sheet "Fig 6D and E" makes paper figure 6D.

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

#isolate cell line total distance mean of each replicate
WT_dist_avg = [WT_rep1.Distance.mean(), WT_rep2.Distance.mean(), WT_rep3.Distance.mean()]
KO_dist_avg = [KO_rep1.Distance.mean(), KO_rep2.Distance.mean(), KO_rep3.Distance.mean()]
KO2_dist_avg = [KO2_rep1.Distance.mean(), KO2_rep2.Distance.mean(), KO2_rep3.Distance.mean()]

#isolate cell line x distance mean of each replicate
WT_x_dist_avg = [WT_rep1['x Distance'].mean(), WT_rep2['x Distance'].mean(), WT_rep3['x Distance'].mean()]
KO_x_dist_avg = [KO_rep1['x Distance'].mean(), KO_rep2['x Distance'].mean(), KO_rep3['x Distance'].mean()]
KO2_x_dist_avg = [KO2_rep1['x Distance'].mean(), KO2_rep2['x Distance'].mean(), KO2_rep3['x Distance'].mean()]

#isolate cell line y distance mean of each replicate
WT_y_dist_avg = [WT_rep1['y Distance'].mean(), WT_rep2['y Distance'].mean(), WT_rep3['y Distance'].mean()]
KO_y_dist_avg = [KO_rep1['y Distance'].mean(), KO_rep2['y Distance'].mean(), KO_rep3['y Distance'].mean()]
KO2_y_dist_avg = [KO2_rep1['y Distance'].mean(), KO2_rep2['y Distance'].mean(), KO2_rep3['y Distance'].mean()]

#combine distance averages
tot_dist_avg = np.concatenate([WT_dist_avg, KO_dist_avg, KO2_dist_avg])
x_dist_avg = np.concatenate([WT_x_dist_avg, KO_x_dist_avg, KO2_x_dist_avg])
y_dist_avg = np.concatenate([WT_y_dist_avg, KO_y_dist_avg, KO2_y_dist_avg])

#make color palette
mypal = {'WT': '#686868', 'KO': '#A4A4A4', 'KO2':'#B9BFBF'}
WTcol = '#686868'
KOcol = '#A4A4A4'
KO2col = '#B9BFBF'

#reformat data for plotting
WT = df[df['Cell Type'] == 'Wild Type']
WTcat = pd.DataFrame()
WTcat['Distance'] = pd.concat([WT['Distance'], WT['x Distance'], WT['y Distance']])
WTcat['Distance Type'] = ['total' for i in range(len(WT['Distance']))] + ['x' for i in range(len(WT['x Distance']))] + ['y' for i in range(len(WT['y Distance']))]
WTcat['Cell Type'] = ['WT' for i in range(len(WTcat))]

KO = df[df['Cell Type'] == 'WASP KO']
KOcat = pd.DataFrame()
KOcat['Distance'] = pd.concat([KO['Distance'], KO['x Distance'], KO['y Distance']])
KOcat['Distance Type'] = ['total' for i in range(len(KO['Distance']))] + ['x' for i in range(len(KO['x Distance']))] + ['y' for i in range(len(KO['y Distance']))]
KOcat['Cell Type'] = ['KO' for i in range(len(KOcat))]

KO2 = df[df['Cell Type'] == 'WASP KO 2']
KO2cat = pd.DataFrame()
KO2cat['Distance'] = pd.concat([KO2['Distance'], KO2['x Distance'], KO2['y Distance']])
KO2cat['Distance Type'] = ['total' for i in range(len(KO2['Distance']))] + ['x' for i in range(len(KO2['x Distance']))] + ['y' for i in range(len(KO2['y Distance']))]
KO2cat['Cell Type'] = ['KO2' for i in range(len(KO2cat))]

dfcat = pd.concat([WTcat, KOcat, KO2cat])

#compare WT, WKO1, WKO2 distances (tot, x, and y) side-by-side
plt.figure()
ax = sns.barplot(x = 'Distance Type', y = 'Distance', hue = 'Cell Type', data = dfcat, palette = mypal)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()

#plot replicate means
plt.scatter(-0.25, tot_dist_avg[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(-0.25, tot_dist_avg[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(-0.25, tot_dist_avg[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0, tot_dist_avg[3], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0, tot_dist_avg[4], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0, tot_dist_avg[5], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.25, tot_dist_avg[6], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.25, tot_dist_avg[7], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.25, tot_dist_avg[8], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.75, x_dist_avg[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.75, x_dist_avg[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.75, x_dist_avg[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(1, x_dist_avg[3], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1, x_dist_avg[4], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1, x_dist_avg[5], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(1.25, x_dist_avg[6], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.25, x_dist_avg[7], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.25, x_dist_avg[8], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(1.75, y_dist_avg[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.75, y_dist_avg[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.75, y_dist_avg[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(2, y_dist_avg[3], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(2, y_dist_avg[4], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(2, y_dist_avg[5], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(2.25, y_dist_avg[6], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(2.25, y_dist_avg[7], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(2.25, y_dist_avg[8], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

#format plot
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


#add statistic brackets for each region (total, x, y)
x1, x2 = -0.25, 0
y, h, col = np.max(tot_dist_avg) + 10, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

x1, x2 = -0.25, 0.25
y, h, col = np.max(tot_dist_avg) + 30, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

#########

x1, x2 = 0.75, 1
y, h, col = np.max(x_dist_avg) + 10, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

x1, x2 = 0.75, 1.25
y, h, col = np.max(x_dist_avg) + 30, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

#######

x1, x2 = 1.75, 2
y, h, col = np.max(y_dist_avg) + 10, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)


x1, x2 = 1.75, 2.25
y, h, col = np.max(y_dist_avg) + 30, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

###statistics###
#paired t-test on replicate means

#WT total distance v KO1/2 total distance
s,p_tot_WT_KO = stats.ttest_rel(WT_dist_avg, KO_dist_avg)
s,p_tot_WT_KO2 = stats.ttest_rel(WT_dist_avg, KO2_dist_avg)
print('p total distance between WT v KO: ' + str(p_tot_WT_KO))
print('p total distance between WT v KO2: ' + str(p_tot_WT_KO2))

#WT x distance v KO1/2 x distance
s,p_x_dist_WT_KO = stats.ttest_rel(WT_x_dist_avg, KO_x_dist_avg)
s,p_x_dist_WT_KO2 = stats.ttest_rel(WT_x_dist_avg, KO2_x_dist_avg)
print('p x distance between WT v KO: ' + str(p_x_dist_WT_KO))
print('p x distance between WT v KO2: ' + str(p_x_dist_WT_KO2))

#WT y distance v KO1/2 y distance
s,p_y_dist_WT_KO = stats.ttest_rel(WT_y_dist_avg, KO_y_dist_avg)
s,p_y_dist_WT_KO2 = stats.ttest_rel(WT_y_dist_avg, KO2_y_dist_avg)
print('p y distance between WT v KO: ' + str(p_y_dist_WT_KO))
print('p y distance between WT v KO2: ' + str(p_y_dist_WT_KO2))


#print mean and standard error
print('Total distance WT:')
print(str(np.round(np.mean(tot_dist_avg[0:3]),1)) + ' ± ' + str(np.round(stats.sem(tot_dist_avg[0:3]),1)))
print('Total distance KO:')
print(str(np.round(np.mean(tot_dist_avg[3:6]),1)) + ' ± ' + str(np.round(stats.sem(tot_dist_avg[3:6]),1)))
print('Total distance KO2:')
print(str(np.round(np.mean(tot_dist_avg[6:]),1)) + ' ± ' + str(np.round(stats.sem(tot_dist_avg[6:]),1)))

print('x distance WT:')
print(str(np.round(np.mean(x_dist_avg[0:3]),1)) + ' ± ' + str(np.round(stats.sem(x_dist_avg[0:3]),1)))
print('x distance KO:')
print(str(np.round(np.mean(x_dist_avg[3:6]),1)) + ' ± ' + str(np.round(stats.sem(x_dist_avg[3:6]),1)))
print('x distance KO2:')
print(str(np.round(np.mean(x_dist_avg[6:]),1)) + ' ± ' + str(np.round(stats.sem(x_dist_avg[6:]),1)))

print('y distance WT:')
print(str(np.round(np.mean(y_dist_avg[0:3]),1)) + ' ± ' + str(np.round(stats.sem(y_dist_avg[0:3]),1)))
print('y distance KO:')
print(str(np.round(np.mean(y_dist_avg[3:6]),1)) + ' ± ' + str(np.round(stats.sem(y_dist_avg[3:6]),1)))
print('y distance KO2:')
print(str(np.round(np.mean(y_dist_avg[6:]),1)) + ' ± ' + str(np.round(stats.sem(y_dist_avg[6:]),1)))
