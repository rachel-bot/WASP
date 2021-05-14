#R. Brunetti, 210331
#Use a paired t-test on replicate means to compare the distributions of
#correlation coefficients between WASP and Arp3 and a negative control where
#the Arp3 channel has been rotated 90 degrees.
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
import statsmodels.stats.api as sms
plt.ion()

# load in correlation coefficients
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 5B')

#separate out replicates
df_1 = df[df.Replicate == 1]
df_2 = df[df.Replicate == 2]
df_3 = df[df.Replicate == 3]

#plot swarmplot
plt.figure(figsize = (3,4))
plt.ylim([-0.2,1])
plt.hlines(0, -1, 2, linewidth = 1, zorder = 1)
ax = sns.swarmplot(x = 'Class', y = 'CorrCoeff', data = df, size=5, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1, alpha = 1)
plt.tight_layout()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#overlay means
plt.scatter(-0.15, df_1[df_1.Class == 'no rotation'].CorrCoeff.mean(), linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0, df_2[df_2.Class == 'no rotation'].CorrCoeff.mean(), linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.15, df_3[df_3.Class == 'no rotation'].CorrCoeff.mean(), linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.85, df_1[df_1.Class == '90 degree rot'].CorrCoeff.mean(), linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1, df_2[df_2.Class == '90 degree rot'].CorrCoeff.mean(), linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.15, df_3[df_3.Class == '90 degree rot'].CorrCoeff.mean(), linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)

#combine means for stats
means_norot = [df_1[df_1.Class == 'no rotation'].CorrCoeff.mean(), df_2[df_2.Class == 'no rotation'].CorrCoeff.mean(), df_3[df_3.Class == 'no rotation'].CorrCoeff.mean()]
means_rot = [df_1[df_1.Class == '90 degree rot'].CorrCoeff.mean(), df_2[df_2.Class == '90 degree rot'].CorrCoeff.mean(), df_3[df_3.Class == '90 degree rot'].CorrCoeff.mean()]

#perform and display stats
statistic, pvalue = scipy.stats.ttest_rel(means_norot, means_rot)
x1, x2 = 0, 1
y, h, col = df['CorrCoeff'].max() + 0.05, 0.05, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*2, "P = "+str(float(round(pvalue, 3))), ha='center', va='bottom', color=col)
print('pvalue: ' + str(pvalue))
