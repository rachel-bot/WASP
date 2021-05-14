# R. Brunetti, 210331
# Script to compare WASP intensities at the beginning of the movie
# when the invagination is open (t_i) with WASP intensities at the time
# of closing (t_f). Data comes from 9 closing invaginations.
# A paired t-test is used to compare intensity at the same invagination
# before and after closing.

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
plt.ion()

# load in WASP intensities and invagination diameters over time
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 3D')

#for each bead put data into a list of arrays
samples = set(df.Bead)
concatdata = []
for sample in samples:
    intensity = df[df.Bead == sample].NormIntensity
    width = df[df.Bead == sample].Diameter
    trace = [int(sample) for i in intensity]
    concatdata.append(np.vstack([trace, width, intensity]))


#identify closing position
#will usually be open --> closed but we built in some safeguards in case there is
#a single frame that was briefly closed or accidentally indicated as closed.
#Currently it will not catch if this occurs more than once, so make sure to cross
#check returned closing frame with FIJI.
closed_pos = []
for trace in concatdata:
    diameter = trace[1]
    indices = [i for i in range(len(diameter)) if diameter[i] == 0] #all closed positions
    offsets = [indices[i]-indices[i-1] for i in range(1, len(indices))] #identify regions that are not consecutive
    if len(offsets)>3: #something >1 to ensure actually closed and not just a single reading of zero (which would likely be error)
        nonzero_ind = indices[np.max(np.argwhere(offsets!=1))+1] #go one up from first closed frame or any closing then opened.
        closed_pos.append(nonzero_ind)
    else:
        closed_pos.append(np.NaN)

#determine intensity at time zero and time of closing (calculated above)
int_f = [concatdata[i][2][closed_pos[i]] for i in range(len(concatdata))]
int_i = [concatdata[i][2][0] for i in range(len(concatdata))]

#convert to dataframe for seaborn plotting
df_int = pd.DataFrame({})
df_int['intensity'] = np.concatenate([int_i, int_f])
df_int['names'] = ['initial' for i in int_i] + ['at closing' for i in int_f]

#plot scatter of pre- and post-closing WASP values and indicate paired measurements
plt.figure()
ax = sns.scatterplot(x = 'names', y = 'intensity', data = df_int, s = 50, color = 'silver', edgecolor = 'k', linewidth = 1, zorder = 10)
for i in range(len(int_i)):
    plt.plot( [0,1], [int_i[i], int_f[i]], c='k', linewidth = 1, zorder = 1)
plt.xlim([-0.2,1.2])
plt.ylim([0,1.15])
ax.set_aspect(3)

#calculate p value from paired t-test on intensities pre- and post- closing
_, p = stats.ttest_rel(int_i, int_f)

# display p value
x1, x2 = 0, 1
y, h, col = df_int.intensity.max()+0.03, 0.04, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h*2-0.03, "P = "+str(float(round(p, 3))), ha='center', va='bottom', color=col)

#print means and standard errors
print('Signal at frame 0: ' + str(np.round(np.mean(int_i),2)) + ' ± ' + str(np.round(stats.sem(int_i),2)))
print('Signal at frame closing: ' + str(np.round(np.mean(int_f),2)) + ' ± ' + str(np.round(stats.sem(int_f),2)))
