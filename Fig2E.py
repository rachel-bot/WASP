# R. Brunetti, 210331
# Read in a csv file of WASP intensities at beads. Compare replicate means between
# total and surface area normalized replicate values using an unpaired t-test.
# There was not a significant difference between total intensities at 100 and 200 nm
# beads despite 200 nm beads having 4x the surface area. When adjusted to account
# for the size difference, intensity per unit area is significantly greater at
# 100 nm beads (3x). Surface area normalized data is presented in Fig 2E.
import os
import skimage
from skimage import io
import numpy as np
import math
import csv
from scipy import optimize,stats
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.stats.api as sms
plt.ion()

# Load data of 2D Gaussian volume of WASP signal at beads
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 2E')

df100 = df[df.Diameter == 100]
df200 = df[df.Diameter == 200]

#how many measurements per experimental replicate?
print('Observations per rep. (100 nm):' + str([len(df100[df100.Replicate == x]) for x in range(1, np.max(df100.Replicate)+1)]))
print('Observations per rep. (200 nm):' + str([len(df200[df200.Replicate == x]) for x in range(1, np.max(df200.Replicate)+1)]))

# separate bead sizes
onehundred_intvals = df100.Intensity
twohundred_intvals = df200.Intensity

# calculate the replicate averages
int_avg =[np.mean(df100[df100.Replicate == 1].Intensity), np.mean(df100[df100.Replicate == 2].Intensity), np.mean(df100[df100.Replicate == 3].Intensity), np.mean(df100[df100.Replicate == 4].Intensity), np.mean(df200[df200.Replicate == 1].Intensity), np.mean(df200[df200.Replicate == 2].Intensity), np.mean(df200[df200.Replicate == 3].Intensity), np.mean(df200[df200.Replicate == 4].Intensity)]
names = ['100', '100', '100', '100', '200', '200', '200', '200']
rep = ['1', '2', '3', '4', '1', '2', '3', '4']
df_avg = pd.DataFrame({'Diameter':names, 'Int':int_avg, 'rep':rep})

#plot raw data
plt.subplot(1,2,1)
plt.title('Total intensities')
ax = sns.swarmplot(x = 'Diameter', y = 'Intensity', color = 'silver', data = df, size=5, edgecolor='gray', linewidth=1, zorder = 1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#perform an unpaired t-test on replicate means
_,p = stats.ttest_ind(df_avg[df_avg.Diameter == '100'].Int, df_avg[df_avg.Diameter == '200'].Int)
print('p-value comparing total intensities:' + str(p))

#overlay replicate means
plt.scatter(-0.15, int_avg[0], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(-0.05, int_avg[1], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.05, int_avg[2], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.15, int_avg[3], linewidth = 1, color = '#FF00FF', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.85, int_avg[4], linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.95, int_avg[5], linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.05, int_avg[6], linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.15, int_avg[7], linewidth = 1, color = '#FF00FF', edgecolor = 'k', zorder = 200, s = 75)

# display p value
x1, x2 = 0, 1
y, h, col = df.Intensity.max(), 5000, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*2-5000, "P = "+str(float(round(p, 3))), ha='center', va='bottom', color=col)

# add field for surface area normalized values to the data frame
df['SA'] = np.concatenate([onehundred_intvals/(4*np.pi*50**2), twohundred_intvals/(4*np.pi*100**2)])
df_avg['SA_val'] = [(4*np.pi*50**2), (4*np.pi*50**2), (4*np.pi*50**2), (4*np.pi*50**2), (4*np.pi*100**2), (4*np.pi*100**2), (4*np.pi*100**2), (4*np.pi*100**2)]
df_avg['SA_norm'] = df_avg['Int']/df_avg['SA_val']

# re-establish these variables given the addition of new fields
df100 = df[df.Diameter == 100]
df200 = df[df.Diameter == 200]

#plot surface area normalized values
plt.subplot(1,2,2)
plt.title('Surface area normalized intensities')
ax = sns.swarmplot(x = 'Diameter', y = 'SA', data = df, color = 'silver', size=5, edgecolor='gray', linewidth=1, zorder = 1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#perform an unpaired t-test on replicate means
_,p_sa = stats.ttest_ind(df_avg[df_avg.Diameter == '100'].SA_norm, df_avg[df_avg.Diameter == '200'].SA_norm)
print('p-value comparing surface area normalized intensities:' + str(p_sa))

#overlay replicate means
plt.scatter(-0.15, df100[df100.Replicate == 1].SA.mean() , linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(-0.05, df100[df100.Replicate == 2].SA.mean() , linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.05, df100[df100.Replicate == 3].SA.mean() , linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.15, df100[df100.Replicate == 4].SA.mean() , linewidth = 1, color = '#FF00FF', edgecolor = 'k', zorder = 200, s = 75)

plt.scatter(0.85, df200[df200.Replicate == 1].SA.mean() , linewidth = 1, color = '#E79D06', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(0.95, df200[df200.Replicate == 2].SA.mean() , linewidth = 1, color = '#53B7E8', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.05, df200[df200.Replicate == 3].SA.mean() , linewidth = 1, color = '#009F72', edgecolor = 'k', zorder = 200, s = 75)
plt.scatter(1.15, df200[df200.Replicate == 4].SA.mean() , linewidth = 1, color = '#FF00FF', edgecolor = 'k', zorder = 200, s = 75)

# display p value
x1, x2 = 0, 1
y, h, col = df100.SA.max() + 0.1, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*2-0.1, "P = "+str(float(round(p_sa, 3))), ha='center', va='bottom', color=col)

#print mean and standard errors
means_100_SANorm = [df100[df100.Replicate == 1].SA.mean(), df100[df100.Replicate == 2].SA.mean(), df100[df100.Replicate == 3].SA.mean(), df100[df100.Replicate == 4].SA.mean()]
means_200_SANorm = [df200[df200.Replicate == 1].SA.mean(), df200[df200.Replicate == 2].SA.mean(), df200[df200.Replicate == 3].SA.mean(), df200[df200.Replicate == 4].SA.mean()]
print('mean intensity of surface area normalized 100nm beads: ' + str(np.round(np.mean(means_100_SANorm),2)) + ' ± ' + str(np.round(stats.sem(means_100_SANorm),2)))
print('mean intensity of surface area normalized 200nm beads: ' + str(np.round(np.mean(means_200_SANorm),2)) + ' ± ' + str(np.round(stats.sem(means_200_SANorm),2)))
