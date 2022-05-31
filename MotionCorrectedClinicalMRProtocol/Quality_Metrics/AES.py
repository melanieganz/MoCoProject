'''
Simon Chemnitz-Thomsen's code to calculate the metric Average Edge Strength. 


Code is based on the article:
Quantitative framework for prospective motion correction evaluation
Nicolas Pannetier, Theano Stavrinos, Peter Ng, Michael Herbst, 
Maxim Zaitsev, Karl Young, Gerald Matson, and Norbert Schuff '''

import numpy as np
from skimage.feature import canny
from scipy.ndimage import convolve
from Quality_Metrics.utils import crop_img, bin_img


def aes(img, brainmask = None, sigma=np.sqrt(2), n_levels = 128, bin = False, crop = True, weigt_avg = False):
    '''
    Parameters
    ----------
    img : numpy array
        Image for which the metrics should be calculated.
    sigma : float
        Standard deviation of the Gaussian filter used 
        during canny edge detection.
    n_levels : int
        Levels of intensities to bin image by
    bin : bool
        Whether or not to bin the image
    crop : bool 
        Whether or not to crop image/ delete empty slices 
        
    Returns
    -------
    AES : float
        Average Edge Strength measure of the input image.
    '''
    #Apply brainmask if given one
    if brainmask is not None: #alternative type(brainmask) != type(None)
        img = img*brainmask
    #Crop image if crop is True
    if crop:
        img = crop_img(img)
    #Bin image if bin is True
    if bin:
        img = bin_img(img, n_levels = n_levels)
    #Centered Gradient kernel in the x-direction
    x_kern = np.array([[-1,-1,-1],
                       [0,0,0],
                       [1,1,1]])
    #Centered Gradient kernel in the y-direction
    y_kern = x_kern.T

    #Shape of volume/img
    vol_shape = np.shape(img)

    #Empty array to contain edge strenghts
    #Function returns the mean of this list
    es = []

    #weights for each slice
    #proportion of non zero pixels
    weights = []

    #Convert to float image
    img = img.astype(float)  # deleted np.

    #For each slice calcule the edge strength
    for slice in range(vol_shape[2]):
        #Slice to do operations on
        im_slice = img[:,:,slice]

        #Weight, proportion of non zero pixels
        weights.append(np.mean(im_slice>0))

        #Convolve slice
        x_conv = convolve(im_slice, x_kern)
        y_conv = convolve(im_slice, y_kern)
        #Canny edge detector
        canny_img = canny(im_slice, sigma = sigma)
        #Numerator and denominator, to be divided
        #defining the edge strength of the slice
        numerator = np.sum(canny_img*( x_conv**2 + y_conv**2 ))
        denominator = np.sum(canny_img)

        #Calculate edge strength
        frac = np.sqrt(numerator)/denominator

        #Append the edge strength
        es.append(frac)
    es = np.array(es)
    #Remove nans
    es  = es[~np.isnan(es)]
    #Return the average edge strength
    if weigt_avg:
        return np.average(es, weights = weights)
    else: return np.mean(es)