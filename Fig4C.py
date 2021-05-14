#R. Brunetti, 210331
# Script to look at change in WASP signal at beads as a cell migrates over them.
# Compares the images from the top and bottom montages in paper figure Fig 4A.
# Beads moving towards the cell rear lose signal while beads at the front gain
# WASP. Note the stark line of transition near the cell middle.
# Cropped image is presented in paper Fig 4C.
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from skimage import io, filters, morphology, measure
plt.ion()

#load image montage
imseq = io.imread('Fig4C_movie.tif')
WASP = imseq[:,0]
beads = imseq[:,1]
final_frame = 25 #corresponds to bottom row in montage in (A)

#identify beads under first and last frame
beadthresh = filters.threshold_otsu(beads)
beads0_mask = beads[0]>beadthresh
beadsf_mask = beads[final_frame]>beadthresh

#generate cell masks for each time point
waspthresh = filters.threshold_otsu(WASP)
WASP0_mask = WASP[0]>0.5*waspthresh
WASPf_mask = WASP[final_frame]>0.5*waspthresh

#identify beads under cell in the final time point
bead_intersection = beadsf_mask*WASPf_mask
otherbeads = beadsf_mask*1-bead_intersection*1
otherbeads = morphology.binary_dilation(otherbeads).astype('float')
otherbeads[otherbeads == 1] = np.NaN
bead_intersection = morphology.binary_dilation(bead_intersection)

#calculate WASP on isolated beads
WASP_signal_0 = (WASP[0]*bead_intersection).astype('float')
WASP_signal_f = (WASP[-1]*bead_intersection).astype('float')

#make a green magenta divergent cmap
cdict = {'green':  ((0.0, 0.0, 0.0),   # no green at 0
                  (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 1.0, 1.0)),  # green at 1

        'red': ((0.0, 1.0, 1.0),   # red at zero (+blue = magenta)
                  (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.0, 0.0)),  # no red at 1

        'blue':  ((0.0, 1.0, 1.0),   # blue at zero (+red = magenta)
                  (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.0, 0.0))   # no blue at 1
       }

# Create the colormap using the dictionary
GnMg = colors.LinearSegmentedColormap('GnMg', cdict)
cmap_choice = plt.cm.get_cmap(GnMg)
cmap_choice.set_bad(color = 'silver') #make beads not under cell (which = NaN) silver

#Create the signal difference map with NaN background beads
diffmap = (WASP_signal_f-WASP_signal_0)+otherbeads
#display
plt.imshow(diffmap, cmap = cmap_choice)
plt.clim([-1000,1000]) #seems to capture the range well, though ultimately just interested in  sign
cbar = plt.colorbar()
plt.axis('off')
cbar.set_ticks([])

#overlay a contour on the cell boundary
contours = measure.find_contours(WASPf_mask, 0.8)
cell_contour = contours[np.argmax([len(x) for x in contours]) ]
plt.plot(cell_contour[:,1], cell_contour[:,0], '--', color = 'k', linewidth = 2)
