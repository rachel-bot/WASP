# R. Brunetti, 210331
# Display the Arp3 intensity at beads under the lamellipod in WT and WASP KO
# backgrounds. Apply statistics to compare values. There was a sufficient
# number of observations to support testing on replicate means. Because there
# was an apparent pattern (where one replicate was notably brighter in both
# conditions), a paired t-test was used. Figure corresponds to Fig 5F.
import os
import skimage
from skimage import io
import numpy as np
import math
import scipy
import csv
from scipy import optimize, stats
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
plt.ion()

#load data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 5F')

mypal = {'WT': '#686868', 'KO': '#B9BFBF'} #color palette
replicate_pal = {1:'#E79D06', 2:'#53B7E8', 3:'#009F72'}

#separate out each replicate
df_1 = df[df.Replicate == 1]
df_2 = df[df.Replicate == 2]
df_3 = df[df.Replicate == 3]


#plot data
plt.figure(figsize = (3,4))
ax = sns.swarmplot(x = 'CellType', y = 'ArpIntensity', data = df, palette = mypal, linewidth=1, zorder = 1)
plt.hlines(0, -0.5,1.5, zorder = 0)
plt.xlim([-0.5,1.5])

#make a data frame of replicate means for seaborn plotting
int_avg = [df_1[df_1['CellType'] == 'WT'].ArpIntensity.mean(), df_2[df_2['CellType'] == 'WT'].ArpIntensity.mean(), df_3[df_3['CellType'] == 'WT'].ArpIntensity.mean(), df_1[df_1['CellType'] == 'KO'].ArpIntensity.mean(), df_2[df_2['CellType'] == 'KO'].ArpIntensity.mean(), df_3[df_3['CellType'] == 'KO'].ArpIntensity.mean()]

rep = ['1', '2', '3', '1', '2', '3']
names = ['WT', 'WT', 'WT', 'KO', 'KO', 'KO']
df_avg = pd.DataFrame({'CellType':names, 'ArpIntensity':int_avg, 'Replicate':rep})


#plot replicate means
plt.scatter(-0.07, int_avg[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0, int_avg[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.07, int_avg[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.93, int_avg[3], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1, int_avg[4], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.07, int_avg[5], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)
plt.ylabel('Background subtracted Arp3 signal')

#Determine and apply the appropriate statistics
#Are there enough observations per replicate to do stats on the replicate means? yes!
WT_obs_per_rep = [len(df_1[df_1.CellType == 'WT']), len(df_2[df_2.CellType == 'WT']), len(df_3[df_3.CellType == 'WT'])]
KO_obs_per_rep = [len(df_1[df_1.CellType == 'KO']), len(df_2[df_2.CellType == 'KO']), len(df_3[df_3.CellType == 'KO'])]
print('WT cells per replicate:' + str(np.asarray(WT_obs_per_rep)))
print('KO cells per replicate:' + str(np.asarray(KO_obs_per_rep)))
#paired t-test (paired chose because there seems to be a pattern in replicate mean values (seemingly signficant difference for one of the replicates maintained in both coniditons))
statistic, pvalue = scipy.stats.ttest_rel(df_avg[df_avg['CellType']=='WT'].ArpIntensity, df_avg[df_avg['CellType']=='KO'].ArpIntensity)
print('pvalue : ' + str(pvalue))

#clean up graph and dispay p-value
plt.tight_layout()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
x1, x2 = 0, 1
y, h, col = df['ArpIntensity'].max()+1000, 500, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*2, "P = "+ str(np.round(pvalue, 3)), ha='center', va='bottom', color=col)

#print means and standard errors
WTmeans = int_avg[:3]
KOmeans = int_avg[3:]
print('mean signal in WT cells ' + str(np.round(np.mean(WTmeans),2)) + ' ± ' + str(np.round(stats.sem(WTmeans),2)))
print('mean signal in KO cells ' + str(np.round(np.mean(KOmeans),2)) + ' ± ' + str(np.round(stats.sem(KOmeans),2)))
