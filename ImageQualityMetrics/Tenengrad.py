
'''
Hannah Eichhorn's code to calculate the metric Tenengrad. 

Code is based on the article: 
Krotkov E. Focusing. Int J Comput Vis. 1988; 1(3):223-237
'''

import numpy as np
from scipy.ndimage import sobel


def TG(img, brainmask=[]):
    '''

    Parameters
    ----------
    img : numpy array
        image for which the metrics should be calculated.
    brainmask_fl : numpy array or list, optional
        If a non-empty numpy array is given, this brainmask will be used to
        mask the images before calculating the metrics. If an empty list is
        given, no mask is applied. The default is [].

    Returns
    -------
    tg : float
        Tenengrad measure of the input image.

    '''
    # image needs to be in floating point numbers in order for gradient to be 
    # correctly calculated
    img = img.astype(np.float64)
    

    # calulate gradients:
    grad_x = sobel(img,axis=0,mode='reflect')
    grad_y = sobel(img,axis=1,mode='reflect')
    grad_z = sobel(img,axis=2,mode='reflect')
    nabla_ab = np.sqrt(grad_x**2+grad_y**2+grad_z**2)
    nabla_abs = nabla_ab.flatten()

    # apply flattened brainmask:
    if len(brainmask) > 0:
        nabla_abs = nabla_abs[brainmask.flatten()>0]

    return np.mean(nabla_abs**2) 
