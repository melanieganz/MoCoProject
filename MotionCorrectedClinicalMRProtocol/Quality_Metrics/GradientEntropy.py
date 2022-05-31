
'''
Hannah Eichhorn's code to calculate the metric Gradient Entropy. 

The code is based on:
McGee K, Manduca A, Felmlee J et al. Image metric-based correction 
(autocorrection) of motion effects: analysis of image metrics. J Magn Reson 
Imaging. 2000; 11(2):174-181

'''

import numpy as np
from scipy.stats import entropy
from scipy.ndimage import sobel


def gradent(img, bm=[]):
    '''
    Parameters
    ----------
    img : numpy array
        image for which the metrics should be calculated.
    bm : numpy array or list, optional
        If a non-empty numpy array is given, this brainmask will be used to
        mask the images before calculating the metrics. If an empty list is
        given, no mask is applied. The default is [].

    Returns
    -------
    ge : float
        Gradient Entropy of the input image.
    '''
    
    # image needs to be in floating point numbers in order for gradient to be 
    # correctly calculated
    img = img.astype(float)

    # calulate gradients:
    grad_x = sobel(img,axis=0,mode='reflect')
    grad_y = sobel(img,axis=1,mode='reflect')
    grad_z = sobel(img,axis=2,mode='reflect')
    nabla_ab = np.sqrt(grad_x**2+grad_y**2+grad_z**2) # maybe needs to be normalized
    
    if len(bm) > 0:
        grad = nabla_ab.flatten()[bm.flatten()!=0]
    else:
        grad = nabla_ab.flatten()
    
    _, counts = np.unique(grad, return_counts=True)
    ge = entropy(counts, base=2)
    
    return ge


