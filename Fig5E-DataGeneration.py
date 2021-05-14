# R. Brunetti, 210331
# User selects beads under the lamellipod in a specified frame. A 2D Gaussian
# is attempted to be fit to WASP signal in the region of the selected beads.
# If the curve fit fails, the fit has an r-squared value less than 0.5, or the
# volume of the Gaussian is extremely high (>100000) or negative a value of zero is assigned-- signifying
# that a Gaussian could not be properly fit. Intensity values will only be used
# as a binary (could a Gaussian be fit or not) in analysis (Fig 5E) since we
# are only looking behavior at a snapshot while the intensity can fluctuate
# (explored in Fig 5F)

#Example WT and KO movies included in Fig5EandF folder

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure, morphology, filters
import scipy
from scipy import optimize
import math
plt.ion()

#user defined values
boxradius = 5 #area to isolate around the bead
FOI = 15 #frame of interest you would like to select beads in
sizestep = 200 #dimension of zoomed inset for easier bead selection

#load image
os.chdir('Fig5EandF')
fname = 'Fig5EandF_WT_ex_cropped.tif'

#order = frame, channel, x, y
imseq = io.imread(fname)
beads = imseq[:,1]
WASP = imseq[:,0]


#display the frame of interest
FOI = FOI - 1
RGB = np.zeros((np.shape(imseq)[2], np.shape(imseq)[3],3))
RGB[:,:,0] = WASP[FOI]/np.max(WASP[FOI])
RGB[:,:,1] = 5*beads[FOI]/np.max(5*beads[FOI]) #scaled so can see both channels well

#select the top left corner of the region to zoom in on
plt.figure(figsize = (12,7))
plt.imshow(RGB)
plt.title('choose left hand corner of zoom')
plt.axis('off')
plt.tight_layout()
plt.draw()

#specify number of beads you would like to analyze in this run
#keep number <= 6 so each puncta and signal kymograph can be evaluated
corner = plt.ginput(1)
plt.imshow(RGB[int(corner[0][1]):int(corner[0][1]+sizestep), int(corner[0][0]):int(corner[0][0]+sizestep),:])
plt.title('How many beads?')
plt.draw()
numbeads = input('Enter number of beads')

#select the beads
plt.draw()
bead_centroid = plt.ginput(int(numbeads))
plt.close()

#extract local WASP signal for each bead
WASP_ROI_tot = []
for i in range(len(bead_centroid)):
    WASP_ROI = WASP[FOI,int(bead_centroid[i][1]+corner[0][1]-boxradius):int(bead_centroid[i][1]+corner[0][1]+boxradius), int(bead_centroid[i][0]+corner[0][0]-boxradius):int(bead_centroid[i][0]+corner[0][0]+boxradius)]
    WASP_ROI_tot.append(WASP_ROI)

#Try to fit a 2D Gaussian
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
bead_ims = WASP_ROI_tot
WASPcorrected = np.zeros((len(WASP_ROI_tot), 2*boxradius, 2*boxradius))
signal = np.zeros((len(bead_ims),1))
print('Goodness of fit values:')
for i in range(len(bead_ims)):
    x = np.linspace(0, 2*boxradius-1, 2*boxradius)
    y = np.linspace(0, 2*boxradius-1, 2*boxradius)
    x, y = np.meshgrid(x, y)
    data_noisy = bead_ims[i] #our data to fit
    initial_guess = (np.max(data_noisy),boxradius,boxradius, 3, 3, 0, 1000) #rough guess from data -- can be improved with trackpy centroid fed in
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
        if gof_rsq > 0.5:
            signal[i] = 2 * math.pi * popt[0] * abs(popt[3]) * abs(popt[4])
        else:
            signal[i] = 0
    except RuntimeError:
        print("Error - curve_fit failed")
        print([i])
        data_fitted = np.zeros((1,2*boxradius*2*boxradius))
        signal[i] = 0
    if signal[i] < 0:
        signal[i] = 0
    if signal[i] > 100000:
        signal[i] = 0
    WASPcorrected[i] = data_fitted.reshape(2*boxradius,2*boxradius)


#display raw WASP signal at bead and the Gaussian fit
maxbeads = 6
plt.figure()
for i in range(len(WASP_ROI_tot)):
    plt.subplot(2,maxbeads,i+1)
    plt.imshow(WASP_ROI_tot[i])
    plt.axis('off')
    plt.subplot(2,maxbeads,i+1+maxbeads)
    plt.imshow(WASPcorrected[i])
    plt.axis('off')
plt.suptitle('top = raw signal and bottom = Gaussian fit')

#print values for transfer to Excel
print('Volume of the Gaussian:')
for value in signal:
    print(int(value[0]))
