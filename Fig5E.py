# R. Brunetti, 210331
# Calculate the fraction of beads under the lamellipod that can recruit any Arp3
# in WT and WASP KO backgrounds.Recruitment was determined by whether a 2D
# Gaussian could be fit to Arp3 signal at a bead.
# An unpaired two-tailed t-test was used to compare the proportion of beads that
# recruited Arp3 in each replicate in WT and WASP KO backgrounds.

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
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 5E')

mypal = {'WT': '#686868', 'KO': '#B9BFBF'}

#separate out each replicate to calculate the number of beads
#that recruit WASP for each background (WT or KO)
df_1 = df[df.Replicate == 1]
tot_num_WT_1 = np.sum(df_1.CellType == 'WT')
nnz_WT_1 = np.sum(df_1[df_1.CellType == 'WT'].ArpIntensity != 0)
tot_num_KO_1 = np.sum(df_1.CellType == 'KO')
nnz_KO_1 = np.sum(df_1[df_1.CellType == 'KO'].ArpIntensity != 0)

df_2 = df[df.Replicate == 2]
tot_num_WT_2 = np.sum(df_2.CellType == 'WT')
nnz_WT_2 = np.sum(df_2[df_2.CellType == 'WT'].ArpIntensity != 0)
tot_num_KO_2 = np.sum(df_2.CellType == 'KO')
nnz_KO_2 = np.sum(df_2[df_2.CellType == 'KO'].ArpIntensity != 0)

df_3 = df[df.Replicate == 3]
tot_num_WT_3 = np.sum(df_3.CellType == 'WT')
nnz_WT_3 = np.sum(df_3[df_3.CellType == 'WT'].ArpIntensity != 0)
tot_num_KO_3 = np.sum(df_3.CellType == 'KO')
nnz_KO_3 = np.sum(df_3[df_3.CellType == 'KO'].ArpIntensity != 0)


#calculate percentage recruiting
percent_WT_1 = nnz_WT_1/tot_num_WT_1
percent_WT_2 = nnz_WT_2/tot_num_WT_2
percent_WT_3 = nnz_WT_3/tot_num_WT_3
percents_WT = [percent_WT_1, percent_WT_2, percent_WT_3]

percent_KO_1 = nnz_KO_1/tot_num_KO_1
percent_KO_2 = nnz_KO_2/tot_num_KO_2
percent_KO_3 = nnz_KO_3/tot_num_KO_3
percents_KO = [percent_KO_1, percent_KO_2, percent_KO_3]

#put into a data frame for seaborn plotting
percentages = np.concatenate([percents_WT, percents_KO])
names = ['WT' for i in range(len(percents_WT))] + ['KO' for i in range(len(percents_KO))]
df_percents = pd.DataFrame({'Percent':percentages, 'Cell Type':names})
plt.figure(figsize = (3,4))
ax = sns.barplot(x = 'Cell Type', y = 'Percent', data = df_percents, palette = mypal)

#plot value for each replicate
plt.scatter(-0.15, np.mean(percent_WT_1), linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0, np.mean(percent_WT_2), linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.15, np.mean(percent_WT_3), linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.85, np.mean(percent_KO_1), linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1, np.mean(percent_KO_2), linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.15, np.mean(percent_KO_3), linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

#apply an unparied t-test to percentage value for each replicate
#There are only three data points, so can't tell if they are normally distributed
#arcsine (or logit) transformation can make percentage data normal
#http://psych.colorado.edu/~carey/Courses/PSYC5741/handouts/Transformations.pdf
transformed_data = np.arcsin(df_percents.Percent)
s, p = stats.ttest_ind(transformed_data[:3], transformed_data[3:])
#Alternatively we find that these stats on non-transformed data yields the same p-value
#while p is slightly different for logit transformation but trend still holds (p = 0.002 here)
print('pvalue:' + str(p))

#clean up plot and add pvalue
plt.tight_layout()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

x1, x2 = 0, 1
y, h, col = df_percents.Percent.max()+0.03, 0.04, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*2-0.03, "P = "+str(float(round(p, 5))), ha='center', va='bottom', color=col)
plt.ylabel('Proportion')

#print mean and standard error
print('Mean percentage Arp3 positive in WT cells: ' + str(np.round(np.mean(percents_WT)*100,2)) + ' ± ' + str(np.round(stats.sem(percents_WT)*100,2)))
print('Mean percentage Arp3 positive in KO cells: ' + str(np.round(np.mean(percents_KO)*100,2)) + ' ± ' + str(np.round(stats.sem(percents_KO)*100,2)))
