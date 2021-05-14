#R. Brunetti, 210331
#Integrated intensity of thresholded WASP signal at beads in 3D STED images
#was calculated and normalized to bead surface area. Values were saved to
#Excel and are displayed and compared here. Image corresponds to paper Fig S4F.

import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os
import seaborn as sns
import numpy as np
plt.ion()

#load data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'S4FFig')

#separate values by bead size
twohundred_ints = df[df.radius == 200].normalizedInt
fivehundred_ints = df[df.radius == 500].normalizedInt

#plot data
sns.swarmplot(x = 'radius', y = 'normalizedInt', data = df, size=5, color = 'silver', edgecolor='gray', linewidth=1, zorder = 1)

#number of observations are lower than ideal
print(str(len(twohundred_ints)) + ' 200 nm beads observed')
print(str(len(fivehundred_ints)) + ' 500 nm beads observed')

#however, data is normally distributed, so still okay to use a t-test
s, p_normal_200 = stats.shapiro(twohundred_ints)
if p_normal_200 > 0.05:
    print('200 nm bead surface area normalized intensities are normally distributed')

s, p_normal_500 = stats.shapiro(fivehundred_ints)
if p_normal_500 > 0.05:
    print('500 nm bead surface area normalized intensities are normally distributed')

#due to limited n, performing a t-test on the full distribution
s, p = stats.ttest_ind(twohundred_ints, fivehundred_ints)
print('pvale = ' + str(p))

#display p value
x1, x2 = 0, 1
y, h, col = df.normalizedInt.max() + 1, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*2-0.5, "P = "+ "{:.2e}".format(p), ha='center', va='bottom', color=col)

#print means and standard errors
print('Mean 200 nm beads: ' + str(np.round(np.mean(twohundred_ints),2)) + ' ± ' + str(np.round(stats.sem(twohundred_ints),2)))
print('Mean 500 nm beads: ' + str(np.round(np.mean(fivehundred_ints),2)) + ' ± ' + str(np.round(stats.sem(fivehundred_ints),2)))
