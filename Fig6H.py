# R. Brunetti, 210331
# Display averaged total MSD traces to compare migration between WT and WASP KO
# cell lines on flat and nanoridged substrates. Corresponds to paper figure 6H.

import pickle as pkl
import os
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import seaborn as sns
import numpy as np
plt.ion()

umperpx = 0.65

#load MSD/persistence data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'Fig 6H and I')

#scale MSD to um/px (twice because value is squared)
df['MSD'] = df['MSD'] * umperpx * umperpx

#separate out different cell backgrounds and conditions
flat_data = df[df.condition == 'flat']
ridged_data = df[df.condition == 'ridged']

#plot these conditions across cell lines
plt.figure(figsize = (8,6))
plt.subplot(1,2,1)
plt.title('Total MSD over 30 min on flat substrates')
sns.lineplot(x = 't', y = 'MSD', hue = 'cell_type', data = flat_data)

plt.subplot(1,2,2)
plt.title('Total MSD over 30 min on ridged substrates')
sns.lineplot(x = 't', y = 'MSD', hue = 'cell_type', data = ridged_data)
plt.tight_layout()
