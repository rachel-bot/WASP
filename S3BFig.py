# R. Brunetti, 210331.
# Compare correlation between overexpressed WASP and endogenous FBP17 in 7.5X7.5
# um ROIs. Data is generated as in Fig5B-DataGeneration.py
# Limited observations per replicate prevented normality assumption, leading us
# to use a paired t-test on all data points to compare WASP-FBP17 correlation
# coefficients to WASP-rotated FBP17 correlation coefficients.

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

#load data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'S3BFig')

#reformat data frame for seaborn plotting
CC = np.concatenate([df.norot, df.rot])
names = ['no rotation' for i in df.norot] + ['rotation' for i in df.rot]
df2 = pd.DataFrame({})
df2['CC'] = CC
df2['names'] = names
df2['rep'] = np.concatenate([df.rep, df.rep])

#separate out replicates to calculate replicate means in each conidtion
rep1 =  df2[df2.rep == 1]
rep2 =  df2[df2.rep == 2]
rep3 =  df2[df2.rep == 3]

norot_means = [rep1[rep1.names == 'no rotation'].mean().CC, rep2[rep2.names == 'no rotation'].mean().CC, rep3[rep3.names == 'no rotation'].mean().CC   ]
rot_means = [rep1[rep1.names == 'rotation'].mean().CC, rep2[rep2.names == 'rotation'].mean().CC, rep3[rep3.names == 'rotation'].mean().CC   ]

#plot data
ax = sns.swarmplot(x = 'names', y = 'CC', data = df2, size=5, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1, alpha = 1)
plt.ylim([-0.3,1])
plt.hlines(0, -1, 2, linewidth = 1, zorder = 1)
plt.tight_layout()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#plot replicate means
plt.scatter(-0.15, norot_means[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0, norot_means[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.15, norot_means[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.85, rot_means[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1, rot_means[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.15, rot_means[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

#check that there are enough observations per replicate to do stats on means
print('number observations per experiment: ' + str([np.sum(rep1.names == 'no rotation'), np.sum(rep2.names == 'no rotation'), np.sum(rep3.names == 'no rotation')]))
# observations per experiment are low but okay if the data is normal
# Let's check:
s,p_rep1 = stats.shapiro(rep1[rep1.names == 'no rotation'].CC)
if p_rep1 > 0.05:
    print('replicate 1 is normal')
else:
    print('replicate 1 is not normal')
s,p_rep2 = stats.shapiro(rep2[rep2.names == 'no rotation'].CC)
if p_rep2 > 0.05:
    print('replicate 2 is normal')
else:
    print('replicate 2 is not normal')
s,p_rep3 = stats.shapiro(rep3[rep3.names == 'no rotation'].CC)
if p_rep3 > 0.05:
    print('replicate 3 is normal')
else:
    print('replicate 3 is not normal')

s,p_rep1 = stats.shapiro(rep1[rep1.names == 'rotation'].CC)
if p_rep1 > 0.05:
    print('replicate 1 (rotated) is normal')
else:
    print('replicate 1 (rotated) is not normal')
s,p_rep2 = stats.shapiro(rep2[rep2.names == 'rotation'].CC)
if p_rep2 > 0.05:
    print('replicate 2 (rotated) is normal')
else:
    print('replicate 2 (rotated) is not normal')
s,p_rep3 = stats.shapiro(rep3[rep3.names == 'rotation'].CC)
if p_rep3 > 0.05:
    print('replicate 3 (rotated) is normal')
else:
    print('replicate 3 (rotated) is not normal')

print('not enough observations to ignore lack of normality')
#They were. However, the rotated values are not coming up as normal...
#Next best course of action is a paired t-test on the whole distribution

#Are full distributions normal?
s,p_norot = stats.shapiro(df2[df2.names == 'rotation'].CC)
if p_rep3 > 0.05:
    print('non-rotated data is normal')
else:
    print('non-rotated data is not normal')
s,p_norot = stats.shapiro(df2[df2.names == 'no rotation'].CC)
if p_rep3 > 0.05:
    print('rotated data is normal')
else:
    print('rotated data is not normal')


#perform paired t-test on all points
_, p = scipy.stats.ttest_rel(df.norot, df.rot)
print('p value from paired t-test on full distribution ' + "{:.2e}".format(p))

#display p value on graph
x1, x2 = 0, 1
y, h, col = df2['CC'].max() + 0.04, 0.04, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*1.5, "P = "+ "{:.2e}".format(p), ha='center', va='bottom', color=col)
plt.ylabel('Correlation Coefficient')
ax.set_aspect(4)

#display means and standard error
print('non-rotated means: ' + str(np.round(np.mean(norot_means),2)) + ' ± ' + str(np.round(stats.sem(norot_means),2)))
print('rotated means: ' + str(np.round(np.mean(rot_means),2)) + ' ± ' + str(np.round(stats.sem(rot_means),2)))
