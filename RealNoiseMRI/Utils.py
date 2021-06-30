'''
Utilitary functions.
'''

import numpy as np
import nibabel as nib
import subprocess
import os
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from scipy.ndimage import sobel


def RegisterNifti(sub, task, subm_dir, intern_dir, out_dir):
    '''
    Uses FreeSurfers registration tool mri_robust_register to transform the
    reconstructed nodding scan (submission) into the space of the ground 
    truth scan (still scan - internal).

    Parameters
    ----------
    sub : str
        Subject ID.
    task : str
        Either 't1' or 't2', specifying which task the submission belongs to.
    subm_dir : str
        Path to the directory with all submisisons from this team.
    intern_dir : str
        Path to the (internal) directory in which the ground truth scans lie.
    out_dir : str
        Path to the directory where the output, i.e. the moved reconstruction,
        should be saved.

    Returns
    -------
    0

    '''
    
    if not os.path.exists(out_dir+sub):
         os.makedirs(out_dir+sub)
    
    # define the relevant filenames:
    movImg = subm_dir+sub+'/Recon_nod_'+task+'.nii'
    targImg = intern_dir+sub+'/Scan_still_'+task+'_off.nii'
    regname = out_dir+sub+'/Recon_nod_'+task+'.lta'
    out = out_dir+sub+'/Recon_nod_'+task+'_mov.nii'
    
    # calculate registration transform which moves movImg into the space of
    # targImg:
    subprocess.run('mri_robust_register --mov '+  movImg + ' --dst ' + targImg + ' --lta ' + regname + ' --satit --iscale', shell=True)
    
    # apply the transform:
    subprocess.run('mri_vol2vol --mov ' + movImg + ' --targ ' + targImg + ' --o ' + out + ' --lta ' + regname, shell=True)
                
    return 0


def Tenengrad(img, brainmask_fl=[]):
    '''
    Calculates the Tenengrad measure as average of gradient magnitudes, as
    described in Pertuz et al. (Pattern Recognition, 2013).

    Parameters
    ----------
    img : numpy array
        image for which the metrics should be calculated.
    brainmask_fl : numpy array or list, optional
        If a non-empty numpy array is given, this brainmask will be used to
        mask the images before calculating the metric. If an empty list is
        given, no mask is applied. The default is [].

    Returns
    -------
    tg : float
        Tenengrad measure of the input image.

    '''
    # image needs to be in floating point numbers in order for gradient to be 
    # correctly calculated
    img = img.astype(np.float)

    # calulate gradients:
    grad_x = sobel(img,axis=0,mode='reflect')
    grad_y = sobel(img,axis=1,mode='reflect')
    grad_z = sobel(img,axis=2,mode='reflect')
    nabla_ab = np.sqrt(grad_x**2+grad_y**2+grad_z**2)
    nabla_abs = nabla_ab.flatten()

    # apply flattened brainmask:
    if len(brainmask_fl) > 0:
        nabla_abs = nabla_abs[brainmask_fl>0]

    return np.mean(nabla_abs**2)   


def Comp_Metrics(image, reference, brainmask, whole_image=False, normal=True):
    '''
    Function that computes the image metrics Structural Similarity, 
    Peak-signal-to-noise ratio and Tenengrad for the given image
    with ref as reference image and brainmask used as mask. The lowest slice
    with non-zero entries is excluded from the analysis to avoid issues 
    from the interpolation during the registration process.

    Parameters
    ----------
    img : str
        Filename of the image for which the metrics should be computed.
    ref : str
        Filename of the reference image.
    brainmask: str
        Filename of the binarized brainmask.
    whole_image : bool, optional
        If whole_image=True, the metrics will be computed on the
        whole images, if False, 0 values of the brainmask are exluded. The
        default is False.
    normal: bool, optional
        If True, the images are normalized before calculating the metrics.
        This option is highly recommended.

    Returns
    -------
    ssim : array
        SSIM values.
    psnr : array
        PSNR values.
    tg : arrayp
        Tenengrad measure values.
    '''
    # load the niftis:
    bm = nib.load(brainmask).get_fdata().astype(np.uint16)
    ref = nib.load(reference).get_fdata().astype(np.uint16)
    img = nib.load(image).get_fdata().astype(float)
    
    # Substitute the lowest non-zero slice from the brainmask to avoid issues
    # that might otherwise arrise from registration-related interpolation:
    for i in range(np.shape(bm)[2]):
        if np.any(bm[:,:,i]):
            n = i
            break
    bm[:,:,n] = np.zeros((np.shape(bm)[0], np.shape(bm)[1]))
    
    # If brainmask should be used, flatten the images and mask them:
    bm_fl = bm.flatten()
    if whole_image == False:
        d_ref = ref.flatten()
        data_ref = d_ref[bm_fl>0]
        d_img = img.flatten()
        data_img = d_img[bm_fl>0]
        print('Values calculated on masked images (0 excluded)')
    else:
        data_ref = ref
        data_img = img
        print('Values calculated on whole images (0 included)')
        
    # Normalize the images:
    if normal == True:
        print('Values calculated on normalized images')
        mean_ref = np.mean(data_ref)
        std_ref = np.std(data_ref)
        data_ref = (data_ref-mean_ref)/std_ref
        
        mean_img = np.mean(data_img)
        std_img = np.std(data_img)
        data_img = (data_img-mean_img)/std_img
        
        mean_img = np.mean(img)
        std_img = np.std(img)
        img_n = (img-mean_img)/std_img
        
    else:
        img_n = img
        
    
    # Calculate PSNR, SSIM and TG:
    peak = np.amax(data_ref)
    psnr = peak_signal_noise_ratio(data_ref, data_img, data_range=peak)
    ssim = structural_similarity(data_ref, data_img, data_range=peak,
                                 gaussian_weights=True)           
    tg = Tenengrad(img_n, brainmask_fl=bm_fl)
        

    return ssim, psnr, tg

