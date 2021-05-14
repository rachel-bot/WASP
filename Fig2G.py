# R. Brunetti, 210331
# Compare the distribution of correlation coefficients from
# WASP with beads to that of known curvature sensors with beads.
# WASP exhibits stronger recruitment to beads than FBP17, CLTA,
# and WAVE complex by an unpaired two-tailed t-test on these
# distributions
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
plt.ion()

# load in correlation coefficients
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 2G')

#combine values into a data frame formatted for seaborn plotting
CC = np.concatenate([df.WASP, df.FBP17, df.CLTA, df.Hem1])
names = ['WASP' for i in range(len(df.WASP))] + ['FBP17' for i in range(len(df.FBP17))] + ['CLTA' for i in range(len(df.CLTA))] + ['WAVE' for i in range(len(df.Hem1))]
df_new = pd.DataFrame({'Protein':names, 'Correlation Coefficient':CC})

#display strip plot of data
ax = sns.stripplot(x = 'Protein', y = 'Correlation Coefficient', size = 7, linewidth = 2, edgecolor = 'k', data = df_new)
plt.hlines(0, -0.5, 4, linewidth = 1, zorder = 1)
plt.xlim([-0.5, 3.5])

#how many observations per marker?
# n is high enough to assume normality
print('num WASP measurements = ' + str(np.sum(np.isfinite(df.WASP))))
print('num FBP17 measurements = ' + str(np.sum(np.isfinite(df.FBP17))))
print('num CLTA measurements = ' + str(np.sum(np.isfinite(df.CLTA))))
print('num WAVE complex measurements = ' + str(np.sum(np.isfinite(df.Hem1))))

#unpaired t-test comparing WASP correlation with beads to
#the correlation of the other proteins of interest with beads
_, p_WASP_FBP17 = scipy.stats.ttest_ind(df.WASP, df.FBP17[~np.isnan(df.FBP17)])
_, p_WASP_CLTA = scipy.stats.ttest_ind(df.WASP, df.CLTA[~np.isnan(df.CLTA)])
_, p_WASP_WAVE = scipy.stats.ttest_ind(df.WASP, df.Hem1)

print('WASP vs FBP17 correlation: p = ' + str(p_WASP_FBP17))
print('WASP vs CLTA correlation: p = ' + str(p_WASP_CLTA))
print('WASP vs WAVE complex correlation: p = ' + str(p_WASP_WAVE))

#print mean and standard error
print('Mean correlaation coefficients ± standard error: ')
print('WASP: ' + str(np.round(np.nanmean(df.WASP),2)) + ' ± ' + str(np.round(stats.sem(df[np.isfinite(df.WASP)].WASP),2)))
print('FBP17: ' + str(np.round(np.nanmean(df.FBP17),2)) + ' ± ' + str(np.round(stats.sem(df[np.isfinite(df.FBP17)].FBP17),2)))
print('CLTA: ' + str(np.round(np.nanmean(df.CLTA),2)) + ' ± ' + str(np.round(stats.sem(df[np.isfinite(df.CLTA)].CLTA),2)))
print('Hem1: ' + str(np.round(np.nanmean(df.Hem1),2)) + ' ± ' + str(np.round(stats.sem(df[np.isfinite(df.Hem1)].Hem1),2)))
