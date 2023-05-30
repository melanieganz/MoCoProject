
'''
Code used for the analysis of the ISMRM 2022 abstract "Evaluating the match 
of image quality metrics with radiological assessment in a dataset with and 
without motion artifacts" 

'''

import nibabel as nib
import numpy as np
import argparse
import sys
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from Tenengrad import TG
from AES import aes
from CoEnt import coent
from ImageEntropy import iment
from GradientEntropy import gradent


def Compute_Metric(filename, metric, brainmask_file=False, ref_file=False, 
                   normal=True):
    '''
    

    Parameters
    ----------
    filename : str
        filename of the nifti image which is supposed to be evaluated.
    metric : str
        which metric to calculate. Please choose between 'SSIM', 'PSNR', 
        'Tenengrad', 'AES', 'GradEntropy', 'ImgEntropy', 'CoEnt' and 'all'. 
        If the option 'all' is chosen, all metrics are returned in a list in 
        the order: SSIM, PSNR, Tenengrad, AES, GradEntropy, ImgEntropy, CoEnt.
    brainmask_file : str, optional.
        filename for the corresponding brainmask. If it is set to False, the
        metric will be calculated on the whole image.Tge default is False.
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
    
    if metric == 'all':
        metrics = ['SSIM', 'PSNR', 'Tenengrad', 'AES', 'GradEntropy', 
                   'ImgEntropy', 'CoEnt']
    else:
        metrics = [metric]
    
    res = []
    for m in metrics:
        if m in ['SSIM', 'PSNR']:
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
            if m == 'PSNR':
                res.append(metric_dict[m](data_ref, data_img, 
                                               data_range=peak))
            else:
                res.append(metric_dict[m](data_ref, data_img, 
                                               data_range=peak, 
                                               gaussian_weights=True))
        
        if m in ['Tenengrad', 'GradEntropy', 'ImgEntropy']:
            if normal == True:
                mean_img = np.mean(img)
                std_img = np.std(img)
                img_n = (img-mean_img)/std_img
            else:
                img_n = img
                
            res.append(metric_dict[m](img_n, brainmask))
            
        if m in ['AES', 'CoEnt']:
            if brainmask_file == False:
                brainmask = None
            res.append(metric_dict[m](img, brainmask))
    
    return res


def main(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None, help="The filename to compute QC metrics over")
    parser.add_argument('-m', '--mask', default=None, help="The filename of the brainmask file")
    parser.add_argument('-c', '--metric', default='all', help="The metric to compute")
    args = parser.parse_args()
    #print(args)
    res = Compute_Metric(args.filename, args.metric, brainmask_file=args.mask)
    print(res[0])
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
