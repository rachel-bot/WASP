# R. Brunetti, 210331
# Calculation script on numpy array of extracted TrackMate trajectories to
# determine distance traveled (in 1D and 2D), the persistence raito, and
# mean squared displacement (in 1D and 2D) across WT and WASP KO cell backgrounds.
# Data analyzed in this example comes from numpy arrays representing one replicate
# day with three fields of view per cell background. Five data frames are
# generated-- df (distance properties),df_dir (persistence), df_msd (total msd),
# df_x_msd (x msd) and df_y_msd (y msd). This is repeated for each replicate and
# each substrate condition (flat or patterned). The whole environment was saved
# using the package "dill" and loaded to generate figure panels.
# For data sharing, df_dir was combined across all replicates in the Excel tab
# "Fig 6C" for flat substrates and "Fig 6D and E" for nanoridged substrates.
# Example trajectories (here the variable WT or KO or KO2) were saved for a
# single replicate and combined for both conditions in Excel tab "Fig 6G".
# Finally, the MSD and persistence data frames were combined across cell lines
# and substrate condition in Excel tab "Fig 6H and I".

# Example numpy data is provided for flat and patterned migration on one day of
# experiments. Specify experiment = 'flat' or not
import skimage
from skimage import io, filters
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os
import seaborn as sns
import pandas as pd
import tifffile
import time
import csv
import random
import dill
plt.ion()

#go to data
prefix = '201102'
os.chdir('Fig6-numpy-examples/')
experiment = 'patterned' #enter flat or patterned/anything else
if experiment == 'flat':
    condition = 'flat-'
else:
    condition = ''

tlength = 80 #length of movie (same across all)

# load in all datasets which have been saved as numpy arrays.
# TrackMate data was converted to these numpys using the function
# array_gen in Fig1G-DataGeneration script
WT1 = np.load(prefix+'-KWC-' + condition + '1.npy', allow_pickle = True)
WT2 = np.load(prefix+'-KWC-' + condition + '2.npy', allow_pickle = True)
WT3 = np.load(prefix+'-KWC-' + condition + '3.npy', allow_pickle = True)
WT = np.concatenate([WT1, WT2, WT3]) #combine datasets
#initialize a replicate value for each cell
rep_WT = ['1' for i in range(len(WT1))] + ['2' for i in range(len(WT2))] + ['3' for i in range(len(WT3))]

KO1 = np.load(prefix+'-WKO1-' + condition + '1.npy', allow_pickle = True)
KO2 = np.load(prefix+'-WKO1-' + condition + '2.npy', allow_pickle = True)
KO3 = np.load(prefix+'-WKO1-' + condition + '3.npy', allow_pickle = True)
KO = np.concatenate([KO1, KO2, KO3])
rep_KO = ['1' for i in range(len(KO1))] + ['2' for i in range(len(KO2))] + ['3' for i in range(len(KO3))]

KO21 = np.load(prefix+'-WKO2-' + condition + '1.npy', allow_pickle = True)
KO22 = np.load(prefix+'-WKO2-' + condition + '2.npy', allow_pickle = True)
KO23 = np.load(prefix+'-WKO2-' + condition + '3.npy', allow_pickle = True)
KO2 = np.concatenate([KO21, KO22, KO23])
rep_KO2 = ['1' for i in range(len(KO21))] + ['2' for i in range(len(KO22))] + ['3' for i in range(len(KO23))]


#keep only cells that are present at least half the movie
tlength_cutoff = int(0.5*tlength)+1
WT_keep_ind = [1 if len(entry) > tlength_cutoff else 0 for entry in WT]
KO_keep_ind = [1 if len(entry) > tlength_cutoff else 0 for entry in KO]
KO2_keep_ind = [1 if len(entry) > tlength_cutoff else 0 for entry in KO2]
rep_WT = [rep_WT[i] for i in range(len(rep_WT)) if WT_keep_ind[i] == 1]
rep_KO = [rep_KO[i] for i in range(len(rep_KO)) if KO_keep_ind[i] == 1]
rep_KO2 = [rep_KO2[i] for i in range(len(rep_KO2)) if KO2_keep_ind[i] == 1]

#isolate the traces from the first half of the movie. Want to compare over
#the same time window
WT = [entry[:tlength_cutoff,:] for entry in WT if len(entry) > tlength_cutoff]
KO = [entry[:tlength_cutoff,:] for entry in KO if len(entry) > tlength_cutoff]
KO2 = [entry[:tlength_cutoff,:] for entry in KO2 if len(entry) > tlength_cutoff]


######################################
# Distance Calculations
######################################

#function to calculate total distance
def absdist_2D(data):
    d_tot = []
    for entry in data:
        xi = entry[:-1, 2]
        xf = entry[1:, 2]
        yi = entry[:-1, 3]
        yf = entry[1:, 3]
        d = np.sum(np.sqrt((xf-xi)**2+(yf-yi)**2))
        d_tot.append(d)
    return d_tot

#function to calculate total distance in specified dir (x or y)
def absdist(data, ind):
    d_tot = []
    for entry in data:
        xi = entry[:-1, ind]
        xf = entry[1:, ind]
        d = np.sum(np.sqrt((xf-xi)**2))
        d_tot.append(d)
    return d_tot

#calculate distance traveled by each cell line in x
WT_x_dist = absdist(WT, 2)
KO_x_dist = absdist(KO, 2)
KO2_x_dist = absdist(KO2, 2)

#calculate distance traveled by each cell line in y
WT_y_dist = absdist(WT, 3)
KO_y_dist = absdist(KO, 3)
KO2_y_dist = absdist(KO2, 3)

#calculate total distance traveled by each cell line
WT_dist = absdist_2D(WT)
KO_dist = absdist_2D(KO)
KO2_dist = absdist_2D(KO2)

#combine distances
dist_tot = np.concatenate([WT_dist, KO_dist, KO2_dist])
dist_x = np.concatenate([WT_x_dist, KO_x_dist, KO2_x_dist])
dist_y = np.concatenate([WT_y_dist, KO_y_dist, KO2_y_dist])

#combine distance data into a data frame and display
names = ['Wild Type' for i in range(len(WT_dist))] + ['WASP KO' for i in range(len(KO_dist))] + ['WASP KO 2' for i in range(len(KO2_dist))]
rep = rep_WT + rep_KO + rep_KO2
df = pd.DataFrame({'Distance':dist_tot, 'Cell Type':names, 'x Distance':dist_x, 'y Distance':dist_y, 'rep':rep})
mypal = {'Wild Type':'#61C371', 'WASP KO':'#96978B', 'WASP KO 2':'#96978B'}

plt.figure(figsize = (10,5))
plt.subplot(1,3,1)
sns.barplot(x = 'Cell Type', y = 'Distance', data = df, palette = mypal)
plt.title('total distance')
plt.ylim([0,500])
ax1 = plt.subplot(1,3,2)
plt.title('x distance')
sns.barplot(x = 'Cell Type', y = 'x Distance', data = df, palette = mypal)
plt.ylim([0,500])
ax2 = plt.subplot(1,3,3)
plt.title('y distance')
sns.barplot(x = 'Cell Type', y = 'y Distance', data = df, palette = mypal)
plt.ylim([0,500])
plt.tight_layout()

#plot trajectories
plt.figure(figsize = (10,6))
plt.suptitle('Trajectories')
plt.xlim([0, 1024])
plt.ylim([0, 1024])
plt.subplot(1,3,1)
plt.title('WT')
for i in range(len(WT)):
    plt.plot(WT[i][:,2], WT[i][:,3], 'b')
plt.xlim([0, 1024])
plt.ylim([0, 1024])
plt.subplot(1,3,2)
plt.title('KO1')
for i in range(len(KO)):
    plt.plot(KO[i][:,2], KO[i][:,3], 'r')
plt.xlim([0, 1024])
plt.ylim([0, 1024])
plt.subplot(1,3,3)
plt.title('KO2')
for i in range(len(KO2)):
    plt.plot(KO2[i][:,2], KO2[i][:,3], 'r')
plt.tight_layout()

######################################
# Persistence Calculations
######################################

#function to calculate the directionality (called persistence ratio in the paper)
#note: because this measurement is between time points, no calculation can be done
#at t = 0
def directionality(data):
    dir_tot = []
    for entry in data:
        dir_tot_row = [np.NaN] #persistence at time zero cannot be calculated
        for i in range(2,len(entry)+1):
            posi_x = entry[0:i-1,2]
            posi_y = entry[0:i-1,3]
            posf_x = entry[1:i,2]
            posf_y = entry[1:i,3]
            displacement = np.sqrt((posf_x[-1]-posi_x[0])**2 + (posf_y[-1]-posi_y[0])**2)
            length_of_track = np.sum(np.sqrt((posf_x-posi_x)**2 + (posf_y-posi_y)**2))
            if length_of_track != 0:
                dir_tot_row.append(displacement/length_of_track)
            else:
                dir_tot_row.append(np.NaN)
        if len(dir_tot_row) < tlength:
            dir_tot_row = dir_tot_row + [np.NaN for j in range(tlength_cutoff-len(dir_tot_row))]
        dir_tot.append(dir_tot_row)
    return dir_tot

#calculate persistence for each cell line and combine across replicates and cell lines for data frame
wtdir = directionality(WT)
kodir = directionality(KO)
ko2dir = directionality(KO2)
WTdir_tot = np.concatenate(wtdir)
KOdir_tot = np.concatenate(kodir)
KO2dir_tot = np.concatenate(ko2dir)
dir_tot = np.concatenate([WTdir_tot, KOdir_tot, KO2dir_tot])
#time is scaled by 0.75 because interval between frame is 45s (0.75 of a minute)
times = 0.75*np.concatenate([np.arange(0,tlength_cutoff) for i in range(len(wtdir)+len(kodir)+len(ko2dir))])
names = ['WT' for i in range(len(WTdir_tot))] + ['KO' for i in range(len(KOdir_tot))] + ['KO2' for i in range(len(KO2dir_tot))]
df_dir = pd.DataFrame({'directionality':dir_tot, 't':times, 'cell_type':names})

#display persistence profiles
plt.figure()
plt.title('Persistence of each cell line')
sns.lineplot(x = 't', y = 'directionality', hue = 'cell_type', data = df_dir)

######################################
# MSD Calculations
######################################

# MSD function for 2D
def MSD(data):
    MSD_tot = []
    for entry in data:
        MSD_row = []
        for i in range(1, len(entry)):
            x_t = entry[:, 2]
            x_0 = x_t[0]*np.ones((1, len(x_t)))
            y_t = entry[:, 3]
            y_0 = y_t[0]*np.ones((1, len(y_t)))
            diff = np.sqrt((x_t-x_0)**2+(y_t-y_0)**2)
            diff_sq = diff**2
            MSD = diff_sq
            MSD_row.append(MSD)
        if len(MSD_row) < tlength:
            MSD_row = MSD_row + [np.NaN for j in range(tlength-len(MSD_row))]
        MSD_tot.append(MSD_row)
    return MSD_tot

# MSD function for 1D (x or y specified)
def MSD_1D(data, axis):
    MSD_tot = []
    for entry in data:
        MSD_row = []
        for i in range(1, len(entry)):
            data_t = entry[:, axis]
            data_0 = data_t[0]*np.ones((1, len(data_t)))
            diff = np.sqrt((data_t-data_0)**2)
            diff_sq = diff**2
            #MSD = np.mean(diff_sq)
            MSD = diff_sq
            MSD_row.append(MSD)
        if len(MSD_row) < tlength:
            MSD_row = MSD_row + [np.NaN for j in range(tlength-len(MSD_row))]
        MSD_tot.append(MSD_row)
    return MSD_tot



#calculate and combine x MSD
wtmsd_x = MSD_1D(WT, 2)
wtmsd_x = [x[0][0] for x in wtmsd_x]
komsd_x = MSD_1D(KO, 2)
komsd_x = [x[0][0] for x in komsd_x]
ko2msd_x = MSD_1D(KO2, 2)
ko2msd_x = [x[0][0] for x in ko2msd_x]

WTMSD_x_tot = np.concatenate(wtmsd_x)
KOMSD_x_tot = np.concatenate(komsd_x)
KO2MSD_x_tot = np.concatenate(ko2msd_x)
MSD_tot_x = np.concatenate([WTMSD_x_tot, KOMSD_x_tot, KO2MSD_x_tot])
times_wt = 0.75*np.concatenate([np.arange(0,tlength) for i in range(len(wtmsd_x)+len(komsd_x)+len(ko2msd_x))])
names = ['WT' for i in range(len(WTMSD_x_tot))] + ['KO' for i in range(len(KOMSD_x_tot))] + ['KO2' for i in range(len(KO2MSD_x_tot))]

#diplay x MSD
plt.figure(figsize = (10,5))
plt.title('x MSD')
plt.subplot(1,3,1)
df_x_msd = pd.DataFrame({'MSD':MSD_tot_x, 't':times, 'cell_type':names})
sns.lineplot(x = 't', y = 'MSD', hue = 'cell_type', data = df_x_msd)
plt.xlim([0,14])
plt.ylim([0,20000])


#repeat process for y
wtmsd_y = MSD_1D(WT, 3)
wtmsd_y = [x[0][0] for x in wtmsd_y]
komsd_y = MSD_1D(KO, 3)
komsd_y = [x[0][0] for x in komsd_y]
ko2msd_y = MSD_1D(KO2, 3)
ko2msd_y = [x[0][0] for x in ko2msd_y]


WTMSD_y_tot = np.concatenate(wtmsd_y)
KOMSD_y_tot = np.concatenate(komsd_y)
KO2MSD_y_tot = np.concatenate(ko2msd_y)
MSD_tot_y = np.concatenate([WTMSD_y_tot, KOMSD_y_tot, KO2MSD_y_tot])
times_wt = 0.75*np.concatenate([np.arange(0,tlength) for i in range(len(wtmsd_y)+len(komsd_y)+len(ko2msd_y))])
names = ['WT' for i in range(len(WTMSD_y_tot))] + ['KO' for i in range(len(KOMSD_y_tot))] + ['KO2' for i in range(len(KO2MSD_y_tot))]

#display
plt.subplot(1,3,2)
plt.title('y MSD')
df_y_msd = pd.DataFrame({'MSD':MSD_tot_y, 't':times, 'cell_type':names})
sns.lineplot(x = 't', y = 'MSD', hue = 'cell_type', data = df_y_msd)
plt.xlim([0,14])
plt.ylim([0,20000])


#repeat process for total
wtmsd = MSD(WT)
wtmsd = [x[0][0] for x in wtmsd]
komsd = MSD(KO)
komsd = [x[0][0] for x in komsd]
ko2msd = MSD(KO2)
ko2msd = [x[0][0] for x in ko2msd]

WTMSD = np.concatenate(wtmsd)
KOMSD = np.concatenate(komsd)
KO2MSD = np.concatenate(ko2msd)
MSD_tot = np.concatenate([WTMSD, KOMSD, KO2MSD])
times_wt = 0.75*np.concatenate([np.arange(0,tlength) for i in range(len(wtmsd)+len(komsd)+len(ko2msd))])
names = ['WT' for i in range(len(WTMSD))] + ['KO' for i in range(len(KOMSD))] + ['KO2' for i in range(len(KO2MSD))]

#display
plt.subplot(1,3,3)
plt.title('total MSD')
df_msd = pd.DataFrame({'MSD':MSD_tot, 't':times, 'cell_type':names})
sns.lineplot(x = 't', y = 'MSD', hue = 'cell_type', data = df_msd)
plt.xlim([0,14])
plt.ylim([0,20000])


#since persistence and MSD data frames are the same size, we can combine them for
#easier raw data reporting
df_dir_msd = pd.DataFrame()
df_dir_msd = df_dir
df_dir['MSD'] = df_msd.MSD
df_dir['xMSD'] = df_x_msd.MSD
df_dir['yMSD'] = df_y_msd.MSD
