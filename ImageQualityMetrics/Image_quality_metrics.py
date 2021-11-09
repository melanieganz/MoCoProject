
'''
Code used for the analysis of the ISMRM 2022 abstract "Evaluating the match 
of image quality metrics with radiological assessment in a dataset with and 
without motion artifacts" 

'''

import nibabel as nib
import numpy as np
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from Tenengrad import TG
from AES import aes
from CoEnt import coent
from ImageEntropy import iment
from GradientEntropy import gradent


def Compute_Metric(filename, brainmask_file, metric, ref_file=False, 
                   normal=True):
    '''
    

    Parameters
    ----------
    filename : str
        filename of the nifti image which is supposed to be evaluated.
    brainmask_file : str
        filename for the corresponding brainmask. If it is set to False, the
        metric will be calculated on the whole image.
    metric : str
        which metric to calculate. Please choose between 'SSIM', 'PSNR', 
        'Tenengrad', 'AES', 'GradEntropy', 'ImgEntropy' and 'CoEnt'.
    ref_file : str, optional
        filename for the reference nifti scan which the image is supposed to 
        be compared to. This is only needed for SSIM and PSNR. The default is 
        False.
    normal : bool, optional
        whether the data should be normalized before metric calculation. The 
        default is True.

    Returns
    -------
    res : float
        value of the metric.

    '''
    
    metric_dict = {'SSIM':structural_similarity, 
                   'PSNR':peak_signal_noise_ratio, 
                   'Tenengrad': TG,
                   'AES': aes,
                   'GradEntropy': gradent,
                   'ImgEntropy': iment,
                   'CoEnt': coent}
    
    img = nib.load(filename).get_fdata().astype(np.uint16)
    if brainmask_file != False:
        brainmask = nib.load(brainmask_file).get_fdata().astype(np.uint16)
    else:
        brainmask = []
    if ref_file != False:
        ref = nib.load(ref_file).get_fdata().astype(np.uint16)
    
    if metric in ['SSIM', 'PSNR']:
        brainmask_fl = brainmask.flatten()
        if brainmask_file != False:
            d_ref = ref.flatten()
            ref_ = d_ref[brainmask_fl>0]
            data_ref = ref_
            dat = img.flatten()
            img_ = dat[brainmask_fl>0]
            data_img = img_
        else:
            data_ref = ref
            data_img = img
        
        if normal == True:
            print('Values calculated on normalized images')
            mean_ref = np.mean(data_ref)
            std_ref = np.std(data_ref)
            data_ref = (data_ref-mean_ref)/std_ref
            
            mean_img = np.mean(data_img)
            std_img = np.std(data_img)
            data_img = (data_img-mean_img)/std_img
        
        peak = np.amax(data_ref)
        if metric == 'PSNR':
            res = metric_dict[metric](data_ref, data_img, data_range=peak)
        else:
            res = metric_dict[metric](data_ref, data_img, data_range=peak,
                                      gaussian_weights=True)
    
    if metric in ['Tenengrad', 'GradEntropy', 'ImgEntropy']:
        if normal == True:
            mean_img = np.mean(img)
            std_img = np.std(img)
            img_n = (img-mean_img)/std_img
        else:
            img_n = img
            
        res = metric_dict[metric](img_n, brainmask)
        
    if metric in ['AES', 'CoEnt']:
        if brainmask_file == False:
            brainmask = None
        res = metric_dict[metric](img, brainmask)
    
    
    return res

