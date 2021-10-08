import numpy as np
import nibabel as nib
import os
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from scipy.ndimage import sobel


def GradientMetrics(img, brainmask_fl=[]):
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


def Comp_Metrics(all_img, ref, brainmask, whole_image=False, compute_gr=True, normal=True):
    '''
    Function that computes the image metrics Peak-signal-to-noise ratio,
    Structural Similarity and Tenengrad for all images
    corresponding to the file names all_img with ref as reference image
    and brainmask used to mask

    Parameters
    ----------
    all_img : list
        contaning filenames of all images for which the metrics should be
        computed.
    ref : numpy array
        reference image.
    brainmask: numpy array
        binarized brainmask loaded from nii file
    whole_image : bool, optional
        If whole_image=True, the metrics will be computed on the
        whole images, if False, 0 values of brainmask are exluded. The
        default is False.
    compute_gr: bool, optional
        If True, tg are computed (takes little more time). If False,
        zeros are returned. The default is True.
    normal: bool, optional
        If True, the images are normalized before calculating the metrics.
        This option is highly recommended.

    Returns
    -------
    psnr : array
        PSNR values.
    ssim : array
        SSIM values.
    tg : arrayp
        Tenengrad measure values.
    names : array
        names of the files (alphabetically sorted).
    '''

    psnr = []
    ssim = []
    tg = []
    brainmask_fl = brainmask.flatten()
    if whole_image == False:
        d_ref = ref.flatten()
        ref_ = d_ref[brainmask_fl>0]
        print('Values calculated on masked images (0 excluded)')
    else:
        print('Values calculated on whole images (0 included)')

    for i in range(len(all_img)):
        print(str(i) + ' done')
        img = nib.load(all_img[i]).get_fdata().astype(np.uint16)
        if whole_image == True:
            data_ref = ref
            data_img = img
        else:
            data_ref = ref_
            dat = img.flatten()
            img_ = dat[brainmask_fl>0]
            data_img = img_

        if normal == True:
            print('Values calculated on normalized images')
            mean_ref = np.mean(data_ref)
            std_ref = np.std(data_ref)
            data_ref = (data_ref-mean_ref)/std_ref
            
            mean_img = np.mean(data_img)
            std_img = np.std(data_img)
            data_img = (data_img-mean_img)/std_img
         
        peak = np.amax(data_ref)
        psnr.append(peak_signal_noise_ratio(data_ref, data_img, data_range=peak))
        ssim.append(structural_similarity(data_ref, data_img, data_range=peak,
                                          gaussian_weights=True))
        #normal=False
        if compute_gr == True:
            if normal == True:
                mean_img = np.mean(img)
                std_img = np.std(img)
                img_n = (img-mean_img)/std_img
            else:
                img_n = img
            teg = GradientMetrics(img_n, brainmask_fl=brainmask_fl)
            tg.append(teg)
        else:
            #ge.append(0)
            tg.append(0)

    for i in range(len(all_img)):
        print(all_img[i])
        print('PSNR: ', psnr[i])
        print('SSIM: ', ssim[i])
        print('Tenengrad: ', tg[i])

    # sort the metrics by name of the files:
    names = []
    for i in range(len(all_img)):
        tmp, filename = os.path.split(all_img[i])

        # search for movement type
        test = filename.find('MOCO_OFF')
        if test > 0:
            descr = 'MOCO_OFF_'
        else:
            test = filename.find('MOCO_ON')
            if test > 0:
                descr = 'MOCO_ON_'
            else:
                test = filename.find('MOCO_RETRO_')
                if test > 0:
                    descr = 'MOCO_RETRO_'
                else:
                    print('ERROR: Filename ' + filename + ' does not contain MOCO description')

        test = filename.find('_NOD_')
        if test > 0:
            descr = descr + 'NOD'
        else:
            test = filename.find('_STILL_')
            if test > 0:
                descr = descr + 'STILL'
            else:
                test = filename.find('_SHAKE_')
                if test > 0:
                    descr = descr + 'SHAKE'
                else:
                    test = filename.find('_SHIFT_')
                    if test > 0:
                        descr = descr + 'SHIFT'                            
                    else:
                        print('EROOR: Filename ' + filename + ' does not contain motion description')

        test = filename.find('_RR_')
        if test > 0:
            descr = descr + '_RR'
            
        test = filename.find('_bad_position_')
        if test > 0:
            descr = descr + '_bad_position'
            
        test = filename.find('_RetroMoCoOff')
        if test > 0:
            descr = descr + '_RetroMoCoOff'
            
        test = filename.find('_RetroMoCoOn')
        if test > 0:
            descr = descr + '_RetroMoCoOn'

        test = filename.find('_REAC_ON_')
        if test > 0:
            descr = descr + '_REAC_ON'
        
        test = filename.find('_REAC_OFF_')
        if test > 0:
            descr = descr + '_REAC_OFF'
        
        names.append(descr)

    names = np.array(names)
    psnr, ssim, tg = np.array(psnr), np.array(ssim), np.array(tg)

    ind = np.argsort(names)
    ssim = ssim[ind][::-1]
    psnr = psnr[ind][::-1]
    tg = tg[ind][::-1]
    names = np.sort(names)[::-1]

    return psnr, ssim, tg, names




