# R. Brunetti, 210331
# Display a single image of a WASP positive cell on beads.
# User identifies and selects beads under the lamellipod.
# A 2D Gaussian is fit to WASP signal in the region of the bead.
# Goodness of fit and the volume of the Gaussian fit are returned.
# Output values are combined across cells/replicates/conditions in Excel.

# Example 100 and 200 nm bead images. There are two WASP-positive beads in each
# example cell that were included in the quantification. Beads were selected based
# on their presence under the lamellipod and maintenance in the TIRF plane,
# under the cell (verified in FIJI). Example cells shows a ~2X increase in signal
# at 200 nm beads compared to 100 nm beads despite 200 nm bead surface area
# being 4X more.

import os
import skimage
from skimage import io
import numpy as np
import math
import scipy
from scipy import optimize
import matplotlib.pyplot as plt
plt.ion()

#user defined. How many pixels around the bead would you like to isolate?
box_radius = 5

#load image data
os.chdir('Fig2E')
fname = 'Fig2E_200_ex.tif'

#image order = frame, channel, x, y
imseq = skimage.io.imread(fname)
bead_chan = imseq[0]
wasp_chan = imseq[1]

#make an RGB image
RGB = np.zeros((np.shape(imseq)[1], np.shape(imseq)[2], 3))
RGB[:,:,0] = 10*np.copy(bead_chan)/np.max(10*bead_chan) #scaled so both channels are comparable
RGB[:,:,1] = np.copy(wasp_chan)/np.max(wasp_chan)
RGB[:,:,2] = np.zeros(np.shape(bead_chan))

#display WASP channel alone and the WASP+bead RGB
plt.figure(figsize = (9,7))
plt.subplot(2,1,1)
plt.title('WASP')
plt.imshow(wasp_chan, cmap = 'Greys')
plt.axis('off')
plt.subplot(2,1,2)
plt.title('WASP + beads')
plt.imshow(RGB)
plt.axis('off')
plt.suptitle('Enter number of beads under the lamellipod')
num_beads = input('How many beads?')

#identify the number of beads specified by clicking them
plt.suptitle('Click the beads in the bottom image')
plt.draw()
coords = plt.ginput(int(num_beads))
coords = [[coords[i][0], coords[i][1]] for i in range(len(coords))]
plt.close()

#get bounding box for each bead
bead_ims = []
max_pixel = []
for i in range(len(coords)):
    new_bead = wasp_chan[int(coords[i][1])-box_radius:int(coords[i][1])+box_radius, int(coords[i][0])-box_radius:int(coords[i][0])+box_radius]
    bead_ims.append(new_bead)
    max_pixel.append(np.max(new_bead))

#from stack exchange, formula for 2D Gaussian
def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                        + c*((y-yo)**2)))
    return g.ravel()


# fit a 2D Gaussian to WASP signal at the beads
# store the signal of this Gaussian (volume under the Gaussian)
# and the 2D fit to compare to raw data
WASPcorrected = np.zeros((len(bead_ims), 2*box_radius, 2*box_radius))
signal = np.zeros((len(bead_ims),1))
print('Goodness of fit values:')
for i in range(len(bead_ims)):
    x = np.linspace(0, 2*box_radius-1, 2*box_radius)
    y = np.linspace(0, 2*box_radius-1, 2*box_radius)
    x, y = np.meshgrid(x, y)
    data_noisy = bead_ims[i] #our data to fit
    initial_guess = (9000,box_radius,box_radius, 2, 2, 0, 1000) #rough guess from data -- can be improved with trackpy centroid fed in
    try:
        popt, pcov = scipy.optimize.curve_fit(twoD_Gaussian, (x, y), data_noisy.ravel(), p0=initial_guess, maxfev = 5000)
        data_fitted = twoD_Gaussian((x, y), *popt)
        # residual sum of squares
        data_noisy_r = data_noisy.ravel()
        ybar = np.mean(data_noisy_r)
        ss = np.sum((data_noisy_r - ybar*np.ones(np.shape(data_noisy_r)))**2)
        gof_sse = np.sum((data_noisy_r-data_fitted)**2)
        gof_rsq = 1 - gof_sse/ss
        print(gof_rsq)
        if gof_rsq > 0.2:
            signal[i] = 2 * math.pi * popt[0] * abs(popt[3]) * abs(popt[4])
        else:
            signal[i] = 0
    except RuntimeError:
        print("Error - curve_fit failed")
        print([i])
        data_fitted = np.zeros((1,2*box_radius*2*box_radius))
        signal[i] = 0
    if signal[i] < 0:
        signal[i] = 0
    WASPcorrected[i] = data_fitted.reshape(2*box_radius,2*box_radius)

#compare raw and fit data
plt.figure()
plt.suptitle('Comparing raw WASP to the 2D fit')
plt.subplot(1,2,1)
plt.title('WASP signal at beads')
plt.imshow(np.concatenate(bead_ims))
plt.axis('off')
plt.subplot(1,2,2)
plt.title('Fitted WASP signal at beads')
plt.imshow(np.concatenate(WASPcorrected))
plt.axis('off')


#print intensity values to transfer to csv
print('WASP signal (volume) from fit:')
for entry in signal:
    print(entry[0])
