#R. Brunetti, 210331
#Script to display xy resolution of 3D STED data, calculated and rendered for the
#paper Fig S4E by G. Kockelkoren. Script is only to show how supporting Excel
#data can yield the histogram in the paper figure.

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

#load data
df = pd.read_excel('paper-data-combined.xlsx', sheet_name = 'S4EFig')


#How to slice up data. Values in nm
approx_minobs = 100
approx_maxobs = 300
binsize = 10
bins = np.linspace(approx_minobs,approx_maxobs,1+(approx_maxobs-approx_minobs)/binsize)

#display histogram
plt.hist(df.confocal_in_nm, bins = bins, facecolor = '#7FAFB3', edgecolor = 'dimgray', label = 'confocal')
plt.hist(df['3DSTED_in_nm'], bins = bins, facecolor = '#8772B1', edgecolor = 'dimgray', label = '3D STED')
plt.legend()
plt.xlabel('xy resolution (nm)')
plt.ylabel('counts')
