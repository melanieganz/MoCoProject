'''
Simon Chemnitz-Thomsen's code to calculate the metric Co-occurence entropy

Code is based on the article:
Quantitative framework for prospective motion correction evaluation
Nicolas Pannetier, Theano Stavrinos, Peter Ng, Michael Herbst, 
Maxim Zaitsev, Karl Young, Gerald Matson, and Norbert Schuff'''

import numpy as np
from skimage.feature.texture import greycomatrix
from utils import bin_img, crop_img


def coent3d(img, brainmask = None, n_levels = 128, bin = True, crop = True, supress_zero = True):
    '''
    Parameters
    ----------
    img : numpy array
        Image for which the metrics should be calculated.
    n_levels : int
        Levels of intensities to bin image by
    bin : bool
        Whether or not to bin the image
    crop : bool 
        Whether or not to crop image/ delete empty slices 
    Returns
    -------
    CoEnt : float
        Co-Occurrence Entropy measure of the input image.
    '''
    #Apply brainmask if given one
    if brainmask is not None: #alternative type(brainmask) != type(None)
        img = img*brainmask
    #Crop image if crop is True
    if crop:
        img = crop_img(img)
    #Bin image if bin is True
    if bin:
        img = bin_img(img, n_levels=n_levels)
    #Scale imgage to have intensity values in [0,255]
    img = 255*(img/np.max(img))
    #Convert image to uint8
    #   as greycomatrix prefers uint8 as input
    img = img.astype(np.uint8)

    #Shape of the image/volume
    vol_shape = np.shape(img)

    #Empty matrix that will be the Co-Occurence matrix
    #Note it is 256x256 as that is the shape of the 
    #output of skimage.feature.greycomatrix
    #even though the image is binned
    co_oc_mat = np.zeros((256,256))


    """
    Note: For 3D encoded images the slice axis does not matter.
    """
    #Generate 2d co-ent matrix for each slice
    for i in range(vol_shape[0]):
        #Temporary co-ent matrix
        tmp_co_oc_mat = greycomatrix(img[i,:,:],
                                 distances = [1],
                                 angles = [0*(np.pi/2),
                                           1*(np.pi/2),
                                           2*(np.pi/2),
                                           3*(np.pi/2)])
        #greycomatrix will generate 4d array
        #The value P[i,j,d,theta] is the number of times 
        #that grey-level j occurs at a distance d and 
        #at an angle theta from grey-level i
        #as we only have one distance we just use 
        #tmp_co_oc_mat[:,:,0,:]
        #As we want the total occurence not split on angles
        #we sum over axis 2.
        tmp_co_oc_mat = np.sum(tmp_co_oc_mat[:,:,0,:], axis = 2)
        #add the occurrences to the co-entropy matrix
        co_oc_mat = co_oc_mat + tmp_co_oc_mat
    
    #Generate 2d co-ent matrix for each slice 
    #   to capture co-occurrence in the direction we sliced before
    for j in range(vol_shape[1]):
        #temporary co-ent matrix
        #note only pi,-pi as angles
        tmp_co_oc_mat = greycomatrix(img[:,j,:],
                                 distances = [1], 
                                 angles = [1*(np.pi/2), 
                                           3*(np.pi/2)])
        #greycomatrix will generate 4d array
        #The value P[i,j,d,theta] is the number of times
        #that grey-level j occurs at a distance d and
        #at an angle theta from grey-level i
        #as we only have one distance we just use
        #tmp_co_oc_mat[:,:,0,:]
        #As we want the total occurence not split on angles
        #we sum over axis 2.
        tmp_co_oc_mat = np.sum(tmp_co_oc_mat[:,:,0,:], axis = 2)
        #add the occurrences to the co-entropy matrix
        co_oc_mat = co_oc_mat + tmp_co_oc_mat
    #Divide by 6 to get average occurance
    co_oc_mat = (1/6)*co_oc_mat
    if supress_zero:
        co_oc_mat[0,0] = 0
    #Normalise
    co_oc_mat = co_oc_mat/np.sum(co_oc_mat)
    #Take log2 to get entropy
    log_matrix = np.log2(co_oc_mat)
    #Return the entropy
    return -np.nansum(co_oc_mat*log_matrix)


def coent2d(img, brainmask = None, n_levels = 128, bin = True, crop = True, supress_zero = True):
    #Apply brainmask if given one
    if brainmask is not None: #alternative type(brainmask) != type(None)
        img = img*brainmask
    #Crop image if crop is True
    if crop:
        img = crop_img(img)
    #Bin image if bin is True
    if bin:
        img = bin_img(img, n_levels=n_levels)
    #Scale imgage to have intensity values in [0,255]
    img = 255*(img/np.max(img))
    #Convert image to uint8
    #   as greycomatrix prefers uint8 as input
    img = img.astype(np.uint8)

    #Shape of the image/volume
    vol_shape = np.shape(img)
    ents = []
    for slice in range(vol_shape[2]):
        tmp_co_oc_mat = greycomatrix(img[:,:,slice],
                                 distances = [1],
                                 angles = [0*(np.pi/2),
                                           1*(np.pi/2),
                                           2*(np.pi/2),
                                           3*(np.pi/2)])
        #greycomatrix will generate 4d array
        #The value P[i,j,d,theta] is the number of times
        #that grey-level j occurs at a distance d and
        #at an angle theta from grey-level i
        #as we only have one distance we just use
        #tmp_co_oc_mat[:,:,0,:]
        #As we want the total occurence not split on angles
        #we sum over axis 2.
        tmp_co_oc_mat = np.sum(tmp_co_oc_mat[:,:,0,:], axis = 2)
        if supress_zero:
            tmp_co_oc_mat[0,0] = 0
        tmp_co_oc_mat =  tmp_co_oc_mat/np.sum(tmp_co_oc_mat)
        log_matrix = np.log2(tmp_co_oc_mat)

        ents.append( -np.nansum(tmp_co_oc_mat*log_matrix) )
        


    return np.nanmean(ents)


def coent(img, brainmask = None, n_levels = 128, bin = True, crop = True, supress_zero = True):
    '''
    Parameters
    ----------
    img : numpy array
        Image for which the metrics should be calculated.
    n_levels : int
        Levels of intensities to bin image by
    bin : bool
        Whether or not to bin the image
    crop : bool 
        Whether or not to crop image/ delete empty slices 
    Returns
    
    -------
    CoEnt : float
        Co-Occurrence Entropy measure of the input image.
    '''

    #Check which function to use:

    #Shape of the volume image
    img_vol = np.shape(img)
    
    #Working under the assumption that the 
    #third axis contains the slices, 
    #eg a 2d encoded image
    #would have shape (256,256,k)
    #where k is the slice number
    #additional assumption: 2d encoded seq does not have more than 100 slices
    #and 3d encoded does not have less than 100 slices.
    if img_vol[2]<100:
        return coent2d(img, brainmask, n_levels, bin, crop, supress_zero)
    else: return coent3d(img, brainmask, n_levels, bin, crop, supress_zero)
