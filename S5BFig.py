# R. Brunetti, 210331
# Compare the lifetimes for puncta originating at the cell front versus those
# originating >= 35% of the cell length (determined empirically from two populations
# emerging in a kernel density estimate). Data is related to that quantified in
# paper figure 1G. We find that one replicate failed to accurately represent puncta
# lifetime, due to reduced temporal sampling, which leads to its exclusion here.
# However, puncta appearance and disappearance were largely unaffected. Posiiton
# of appearance, disappearance, and the number of frames a puncta existed were
# all manually calculated using FIJI and TrackMate.
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv
from skimage import io
from scipy import ndimage, stats, optimize
import math
plt.ion()


#load data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'S5BFig')

#separate by means
df1 = df[df.replicate == 1]
df2 = df[df.replicate == 2]
df3 = df[df.replicate == 3]

#Are there enough puncta per cell to do stats on means?-- yes
puncta_per_cell_1 = [np.sum(df1.cell == x) for x in np.unique(df1.cell)]
puncta_per_cell_2 = [np.sum(df2.cell == x) for x in np.unique(df2.cell)]
puncta_per_cell_3 = [np.sum(df3.cell == x) for x in np.unique(df3.cell)]
print(str(puncta_per_cell_1) + ' observations per cell in replicate 1')
print(str(puncta_per_cell_2) + ' observations per cell in replicate 2')
print(str(puncta_per_cell_3) + ' observations per cell in replicate 3')

#calculate cell means for appearance
app_means_1 = [np.mean(df1[df1.cell == x].appearancePos) for x in np.unique(df1.cell)]
app_means_2 = [np.mean(df2[df2.cell == x].appearancePos) for x in np.unique(df2.cell)]
app_means_3 = [np.mean(df3[df3.cell == x].appearancePos) for x in np.unique(df3.cell)]

#calculate cell means for disappearance
disapp_means_1 = [np.mean(df1[df1.cell == x].disappearancePos) for x in np.unique(df1.cell)]
disapp_means_2 = [np.mean(df2[df2.cell == x].disappearancePos) for x in np.unique(df2.cell)]
disapp_means_3 = [np.mean(df3[df3.cell == x].disappearancePos) for x in np.unique(df3.cell)]

#calculate cell means for lifetime
lifetime_means_1 = [np.mean(df1[df1.cell == x].lifetime) for x in np.unique(df1.cell)]
lifetime_means_2 = [np.mean(df2[df2.cell == x].lifetime) for x in np.unique(df2.cell)]
lifetime_means_3 = [np.mean(df3[df3.cell == x].lifetime) for x in np.unique(df3.cell)]

#combine data
rep = [1 for i in app_means_1] + [2 for i in app_means_2] + [3 for i in app_means_3]
appmeans = np.concatenate([app_means_1, app_means_2, app_means_3])
disappmeans = np.concatenate([disapp_means_1, disapp_means_2, disapp_means_3])
lifetimemeans = np.concatenate([lifetime_means_1, lifetime_means_2, lifetime_means_3])

#plot cell means
plt.figure(figsize = (15,4))
plt.subplot(1,3,1)
plt.ylabel('Position of appearance')
plt.xlabel('replicate')
sns.swarmplot(x = rep, y = appmeans)
plt.title('')
plt.ylim([0,1])
plt.subplot(1,3,2)
plt.ylabel('Position of disappearance')
plt.xlabel('replicate')
sns.swarmplot(x = rep, y = disappmeans)
plt.ylim([0,1])
plt.subplot(1,3,3)
plt.ylabel('Normalized lifetime')
plt.xlabel('replicate')
sns.swarmplot(x = rep, y = lifetimemeans/np.max(lifetimemeans))
plt.ylim([0,1])
plt.tight_layout()
plt.suptitle('comparing replicates for goodness of tracking')

#print statistics on cell means
print('p appearance 1 v 2: =' + str(stats.ttest_ind(app_means_1, app_means_2)[1]))
print('p appearance 1 v 3: =' + str(stats.ttest_ind(app_means_1, app_means_3)[1]))
print('p appearance 2 v 3: =' + str(stats.ttest_ind(app_means_2, app_means_3)[1]))

print('p disappearance 1 v 2: =' + str(stats.ttest_ind(disapp_means_1, disapp_means_2)[1]))
print('p disappearance 1 v 3: =' + str(stats.ttest_ind(disapp_means_1, disapp_means_3)[1]))
print('p disappearance 2 v 3: =' + str(stats.ttest_ind(disapp_means_2, disapp_means_3)[1]))

print('p lifetimeearance 1 v 2: =' + str(stats.ttest_ind(lifetime_means_1, lifetime_means_2)[1]))
print('p lifetimeearance 1 v 3: =' + str(stats.ttest_ind(lifetime_means_1, lifetime_means_3)[1]))
print('p lifetimeearance 2 v 3: =' + str(stats.ttest_ind(lifetime_means_2, lifetime_means_3)[1]))

#decision to omit replicate 3, which was used decreased temporal sampling
print('While rep 3 performed similarly to the other replicates for puncta appearance')
print('and disappearance, it poorly recapitulated lifetime and could not be used')
print('in this analysis. Due to decreased temporal sampling in this replicate.')

#combine data from replicate 1 and replicate 2
lifetime = np.concatenate([df1.lifetime, df2.lifetime])
x_dis = np.concatenate([df1.disappearancePos, df2.disappearancePos])
x_app = np.concatenate([df1.appearancePos, df2.appearancePos])
names = ['1' for i in range(len(df1))] + ['2' for i in range(len(df2))]

#display normalized lifetime as a function the appearance position via a kernel density estimate
plt.figure()
ax = sns.kdeplot(x_app, lifetime/np.max(lifetime), shade = True, cmap = 'magma', shade_lowest = False)
plt.xlabel('Relative position of appearance relative to cell front')
plt.ylabel('Normalized lifetime')
plt.xlim([0, 1])
plt.ylim([0, 1])

#make a data frame for plotting
df_rep1and2 = pd.DataFrame()
df_rep1and2['appearance'] = x_app
df_rep1and2['cell'] = np.concatenate([df1.cell, df2.cell])
df_rep1and2['names'] = names
df_rep1and2['lifetime'] = lifetime
#classify puncta as belonging to the front nucleating or back nucleating populations
#that arose in KDE plot
#value of 0.35 determine from KDE plot
df_rep1and2['pos'] = ['front' if x <0.35 else 'back' for x in df_rep1and2.appearance]

#display these two populations in a histogram
plt.figure()
x = plt.hist(df_rep1and2[df_rep1and2.pos == 'front'].lifetime, density = True, edgecolor='k', label = 'front')
plt.hist(df_rep1and2[df_rep1and2.pos == 'back'].lifetime, bins = x[1], density = True, alpha = 0.5, edgecolor='k', label = 'rear')
plt.legend()

#what are the correct stats to do?
#how many front and how many rear observations per cell?
df_rep1 = df_rep1and2[df_rep1and2.names == '1']
front_obs_per_cell_1 = [np.sum(df_rep1[df_rep1.cell == x].pos == 'front') for x in np.unique(df1.cell)]
back_obs_per_cell_1 = [np.sum(df_rep1[df_rep1.cell == x].pos == 'back') for x in np.unique(df1.cell)]
df_rep2 = df_rep1and2[df_rep1and2.names == '2']
front_obs_per_cell_2 = [np.sum(df_rep2[df_rep2.cell == x].pos == 'front') for x in np.unique(df2.cell)]
back_obs_per_cell_2 = [np.sum(df_rep2[df_rep2.cell == x].pos == 'back') for x in np.unique(df2.cell)]
print('front observations per cell: ' + str(np.concatenate([front_obs_per_cell_1, front_obs_per_cell_2])))
print('rear observations per cell: ' + str(np.concatenate([back_obs_per_cell_1, back_obs_per_cell_2])))
#not enough obervations in each region per cell

#So, are there enough observations on the replicate level
front_obs_per_rep_1 = np.sum(df_rep1.pos == 'front')
back_obs_per_rep_1 = np.sum(df_rep1.pos == 'back')
front_obs_per_rep_2 = np.sum(df_rep2.pos == 'front')
back_obs_per_rep_2 = np.sum(df_rep2.pos == 'back')
print('front observations per replicate: ' + str([front_obs_per_rep_1, front_obs_per_rep_2]))
print('back observations per replicate: ' + str([back_obs_per_rep_1, back_obs_per_rep_2]))
#still not enough. Back nucleating observations are sufficiently rare it appears!

#Need to apply stats to full distributions
#sample size supports normality assumption and therefore use of a t-test
s,p = stats.ttest_ind(df_rep1and2[df_rep1and2.pos == 'front'].lifetime, df_rep1and2[df_rep1and2.pos == 'back'].lifetime)
print('p = ' + str(np.round(p,4)))

#print mean and standard error
lifetime_front = np.concatenate([df_rep1[df_rep1.pos == 'front'].lifetime, df_rep2[df_rep2.pos == 'front'].lifetime])
lifetime_back = np.concatenate([df_rep1[df_rep1.pos == 'back'].lifetime, df_rep2[df_rep2.pos == 'back'].lifetime])
print('mean front lifetime: ' + str(0.5*np.round(np.mean(lifetime_front),2)) + ' ± ' + str(0.5*np.round(stats.sem(lifetime_front),2)))
print('mean back lifetime: ' + str(0.5*np.round(np.mean(lifetime_back),2)) + ' ± ' + str(0.5*np.round(stats.sem(lifetime_back),2)))
