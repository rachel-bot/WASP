# R. Brunetti, 210331
# Script to read in and plot WASP distribution percentages.
# Right hand figure is Fig. 1C
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import scipy
import statsmodels.stats.api as sms
plt.ion()

#load data
dataDoG = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 1C')

#make a new, reformatted data frame for seaborn plotting
df = pd.DataFrame({})
df['Region'] = ['Leading Edge' for i in dataDoG.LE] + ['Cell Front' for i in dataDoG.front] + ['Cell Rear' for i in dataDoG.rear] + ['Uropod' for i in dataDoG.U]
df['Replicate'] = np.concatenate([dataDoG.replicate, dataDoG.replicate, dataDoG.replicate, dataDoG.replicate])
df['Percent'] = 100*np.concatenate([dataDoG.LE, dataDoG.front, dataDoG.rear, dataDoG.U])
df['AreaNormalizedPercent'] = 100*np.concatenate([dataDoG.LE_area_norm, dataDoG.front_area_norm, dataDoG.rear_area_norm, dataDoG.U_area_norm])

#plot the raw percent values
plt.figure(figsize = (10,4))
ax1 = plt.subplot(1,2,1)
sns.barplot(x = 'Region', y = 'Percent', data = df)
plt.ylim([0,80])
plt.ylabel('Percent')
ax1.invert_xaxis()

#plot the area-normalized percent values
ax2 = plt.subplot(1,2,2)
sns.barplot(x = 'Region', y = 'AreaNormalizedPercent', data = df)
plt.ylabel('Area Normalized Percent')
plt.ylim([0,80])
ax2.invert_xaxis()
plt.tight_layout()

#separate out regions
LE =  df[df.Region == 'Leading Edge'].AreaNormalizedPercent
front =  df[df.Region == 'Cell Front'].AreaNormalizedPercent
rear =  df[df.Region == 'Cell Rear'].AreaNormalizedPercent
U =  df[df.Region == 'Uropod'].AreaNormalizedPercent

#print area-normalized means and standard error
print('means ± standard error:')
print('Leading Edge: ' + str(np.round(np.mean(LE), 1)) + ' ± ' + str(np.round(stats.sem(LE),1)))
print('Cell Front: ' + str(np.round(np.mean(front), 1)) + ' ± ' + str(np.round(stats.sem(front),1)))
print('Cell Rear: ' + str(np.round(np.mean(rear), 1)) + ' ± ' + str(np.round(stats.sem(rear),1)))
print('Uropod: ' + str(np.round(np.mean(U), 1)) + ' ± ' + str(np.round(stats.sem(U),1)))
