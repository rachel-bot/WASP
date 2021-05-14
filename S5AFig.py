# R. Brunetti, 210331
# Generate single bead traces of normalized intensity as a function of distance
# from the leading edge. Calculated for beads that traverse the majority of the
# cell length. Data is averaged into a discretized line plot (to the
# nearest tenth of a cell) with raw values from the traces overlaid on top and
# colored by the KDE of the intensity distribution in that position bin.
# Creates paper Fig S5A.
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure, morphology, filters
import pandas as pd
import seaborn as sns
import math
from scipy import stats, ndimage
import random
plt.ion()

# load in distance and relative intensity values
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 4B and S5AFig')

#plot all single bead traces
#isolate by replicate
df_1 = df[df.replicate == 1]
df_2 = df[df.replicate == 2]
df_3 = df[df.replicate == 3]

#function to plot all single cell traces in a replicate
def trace_plotter(df_rep):
    rep_cells = np.unique(df_rep.cell)
    for cell in rep_cells:
        df_cell_temp = df_rep[df_rep.cell == cell]
        beads_per_cell = np.unique(df_cell_temp.bead)
        for bead in beads_per_cell:
            bead_data = df_cell_temp[df_cell_temp.bead == bead]
            plt.plot(bead_data.relDist, bead_data.relInt, '-')

#display all single bead traces
plt.figure(figsize = (11,6))
plt.subplot(1,2,1)
trace_plotter(df_1)
trace_plotter(df_2)
trace_plotter(df_3)
plt.xlabel('Relative distance from cell front')
plt.ylabel('Normalized intensity')

#function to discretize traces for grouping
def round_nearest(x, a):
    return np.round(np.round(x / a) * a, -int(math.floor(math.log10(a))))

#discretize position for KDE stripplot
df['Position'] = round_nearest(df.relDist,0.1)

#calculate KDE
x = df.relDist
y = df.relInt
xy = np.vstack([x,y])
z = stats.gaussian_kde(xy)(xy)

#make a density for each positional bin
bins = np.linspace(0,1,11)
for value in bins:
    value = np.round(value,2)
    kernel = stats.gaussian_kde(df[df.Position == value].relInt)
    df.loc[df.Position == value, 'colors'] = kernel(df[df.Position == value].relInt)


#plot x discretized data lineplot and the kde-encoded data distribution
df['d10'] = df.Position*10 #needed to scale to plot lineplot and stripplot on same axes
ax = plt.subplot(1,2,2)
sns.stripplot(x = 'Position', y = 'relInt', hue = 'colors', data = df, linewidth = 1, alpha = 0.7, size = 8, palette = 'magma')
sns.lineplot(x = 'd10', y = 'relInt', data = df, color = 'black', ci = 95, zorder = 10, err_style = 'band', linewidth = 2)
legend = ax.legend()
legend.remove()
plt.xlabel('Relative distance from cell front')
plt.ylabel('Normalized intensity')
plt.tight_layout()
