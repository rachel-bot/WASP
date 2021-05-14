# R. Brunetti, 210331
# Script that combines TrackMate csv of puncta tracking data with
# list of user-determined indices where the time and position of puncta
# appearance was accurately identified. Based off of script for Fig 1G (EDTA experiment).
# Each frame where a puncta appearance occurs is displayed with puncta
# that just appeared marked. User selects cell rear and cell front. A line is drawn
# between these points and the puncta coordinate is projected along it.
# This value is normalized to cell length to get a relative position of
# appearance relative to the cell front. Values are printed for concatenation
# with other cell/replicate values in excel.
# Note : movies are only 10s long so measurements will not be biased by
# consistent reorientation in Cdc42 KO cells. We are capturing short bursts of
# persistent WT-like movement.

# WT example data is cell 10 in replicate 1, which corresponds to rows 2-12.
# Cdc42 KO example data is cell 2 in replicate 1, which corresponds to rows 93-101.
# Values will be approximate as they rely on user user selected positions.

import os
import skimage
from skimage import io, filters, measure, morphology
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import numpy as np
import pickle as pkl
import pandas as pd
import csv
plt.ion()

# function that converts trackmate csv to a list of arrays
# each array is [puncta index, frame num, x coord, y coord]
def array_gen(fname):
    i = 0
    array_rep = []
    with open(fname.split('.')[0] + '.csv', newline='') as csvfile:
         spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
         for row in spamreader:
            i+= 1
            if i ==1:
                title_row = row
            else:
                particle = row[0].split(',')[2]
                x = float(row[0].split(',')[4])
                y = float(row[0].split(',')[5])
                t = row[0].split(',')[8]
                array_rep.append([int(particle),int(t),float(x),float(y)])
    array_rep = np.asarray(array_rep)
    particlenums = np.unique(array_rep[:,0])
    array_rep_total = []
    for particle in particlenums:
        array_rep_total.append(array_rep[array_rep[:,0] == particle])
    return array_rep_total

# go to example data
os.chdir('Fig4F')
nameindiv = 'WT10' #name of cell
experiment = 'WT' #WT or KO
tracks_ex = 'Fig4F_' + nameindiv + '-TM.csv'
indices_ex = 'Fig4F_' + experiment + '-indices-persistent.csv'
movie_ex = 'Fig4F_' + experiment + '_example.tif'

# apply function to condense TrackMate csv to relevant fields
tracked_particles = array_gen(tracks_ex)
# load in cell image sequence
imseq = io.imread(movie_ex)

# Read in master csv for the replicate+condition that has the puncta indices
# where the appearance position was accurately tracked.
# Isolate indices for the cell of interest
indices = pd.read_csv(indices_ex)[nameindiv]
# Keep puncta whose index is in the list of correctly tracked puncta
puncta_kept = [entry for entry in tracked_particles if entry[0][0] in np.asarray(indices)]

# Get [time, x, y] for the first frame of every puncta
# and sort by frame number
firstframe = np.asarray([[x[0][1], x[0][2], x[0][3]] for x in puncta_kept])
firstframe_sorted = firstframe[firstframe[:,0].argsort()]
# Isolate and order frames where a puncta appeared
uniqueframes = set(firstframe_sorted[:,0])
uniqueframes = np.sort([x for x in uniqueframes])

# For each time point where a puncta appears plot a marker on the puncta.
# Manually select the rear and then front of each cell (used to calculate relative position).
# Relative position is calculated by projecting punta coord onto the line from back to front
# and normalizing to the total cell length
distance_along_axis = []
for val in uniqueframes:
    coords = firstframe_sorted[firstframe_sorted[:,0] == val]
    plt.imshow(imseq[int(val)])
    plt.title('Select the cell rear and then cell front')
    for entry in coords:
        plt.plot(entry[1], entry[2], 'r.')
    AB = plt.ginput(2)
    for pt in AB:
        plt.plot(pt[0], pt[1], 'x')
    plt.show()
    v1 = [AB[0][0]- AB[1][0], AB[0][1]-AB[1][1]]


    for entry in coords:
        v2 = [entry[1]- AB[1][0], entry[2]-AB[1][1]]
        distance_along_axis.append(np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v1))

    plt.close()

#plot histogram for visualization
plt.title('Relative position of appearance relative to the cell front')
plt.hist(distance_along_axis, bins = np.linspace(0,1,11)) #discretized to nearest 10% cell length
plt.xlim([0,1])

#print data for transfer to a csv
print('Position of appearance:')
for entry in distance_along_axis:
    print(entry)
