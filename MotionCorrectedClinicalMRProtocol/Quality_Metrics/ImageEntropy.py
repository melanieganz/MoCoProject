
'''
Hannah Eichhorn's code to calculate the metric Image Entropy. 

The code is based on:
McGee K, Manduca A, Felmlee J et al. Image metric-based correction 
(autocorrection) of motion effects: analysis of image metrics. J Magn Reson 
Imaging. 2000; 11(2):174-181
'''

import numpy as np
from scipy.stats import entropy


def iment(img, bm=[]):
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
    ie : float
        Image Entropy of the input image.
    '''
    
    
    # flatten and brainmask the image
    if len(bm) > 0:
        image = img.flatten()[bm.flatten!=0]
    else:
        image = img.flatten()
    
    _, counts = np.unique(image, return_counts=True)
    ie = entropy(counts, base=2)
    
    return ie


