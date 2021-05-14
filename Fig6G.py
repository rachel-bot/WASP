# R. Brunetti, 210331
# Display trajectories of WT and WASP KO cells migrating for one hour on
# flat and nanoridged substrates. Displays data in paper Fig 6G.

import pickle as pkl
import os
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import seaborn as sns
import numpy as np
plt.ion()

#px per image for scaling
image_dim = 1024

#load trajectory data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 6G')

#separate out different cell backgrounds and conditions
WTdata = df[df['cell type']== 'WT']
WT_flat = WTdata[WTdata.condition == 'flat']
WT_ridged = WTdata[WTdata.condition == 'ridged']


KOdata = df[df['cell type']== 'KO']
KO_flat = KOdata[KOdata.condition == 'flat']
KO_ridged = KOdata[KOdata.condition == 'ridged']

#find unique cells in each experiment to plot trajectory
WT_flat_unique_cells = np.unique(WT_flat.cell)
plt.subplot(2,2,1)
plt.xlim([0,1024])
plt.ylim([0,1024])
plt.xticks([])
plt.yticks([])
plt.title('1 hr of WT on flat substrate')
for cellnum in WT_flat_unique_cells:
    plt.plot(WT_flat[WT_flat.cell == cellnum].x, WT_flat[WT_flat.cell == cellnum].y, color = 'purple')

KO_flat_unique_cells = np.unique(KO_flat.cell)
plt.subplot(2,2,2)
plt.xlim([0,1024])
plt.ylim([0,1024])
plt.xticks([])
plt.yticks([])
plt.title('1 hr of KO on flat substrate')
for cellnum in KO_flat_unique_cells:
    plt.plot(KO_flat[KO_flat.cell == cellnum].x, KO_flat[KO_flat.cell == cellnum].y, color = 'skyblue')

WT_ridged_unique_cells = np.unique(WT_ridged.cell)
plt.subplot(2,2,3)
plt.xlim([0,1024])
plt.ylim([0,1024])
plt.xticks([])
plt.yticks([])
plt.title('1 hr of WT on ridged substrate')
for cellnum in WT_ridged_unique_cells:
    plt.plot(WT_ridged[WT_ridged.cell == cellnum].x, WT_ridged[WT_ridged.cell == cellnum].y, color = 'purple')

KO_ridged_unique_cells = np.unique(KO_ridged.cell)
plt.subplot(2,2,4)
plt.xlim([0,1024])
plt.ylim([0,1024])
plt.xticks([])
plt.yticks([])
plt.title('1 hr of KO on ridged substrate')
for cellnum in KO_ridged_unique_cells:
    plt.plot(KO_ridged[KO_ridged.cell == cellnum].x, KO_ridged[KO_ridged.cell == cellnum].y, color = 'skyblue')
