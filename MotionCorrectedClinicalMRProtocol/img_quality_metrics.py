import numpy as np
import nibabel as nib
import os
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from Quality_Metrics.Tenengrad import TG
from Quality_Metrics.AES import aes
from Quality_Metrics.CoEnt import coent
from Quality_Metrics.ImageEntropy import iment
from Quality_Metrics.GradientEntropy import gradent


def Calc_Metrics(filename, metric, brainmask_file=False, ref_file=False, 
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


def Comp_All_Metrics(all_img, bm_file, ref_file, out_dir, save, tag, normal=True):
    '''
    Computes the metrics SSIM, PSNR, Tenengrad, AES, GradEntropy, ImgEntropy, 
    CoEnt for all nifti files specified in all_img. 

    Parameters
    ----------
    all_img : list or array of str.
        list of nifti images.
    bm_file : str
        filename and path for the corresponding brainmask.
    ref_file : str
        filename and path for the corresponding reference image.
    out_dir : str
        path to the output directory where the metrics should be stored.
    save : str
        date tag for naming the output file.
    tag : str
        sequence tag for naming the output file.
    normal : bool, optional
        Whether the images should be normalized. The default is True.

    Returns
    -------
    int
        DESCRIPTION.

    '''
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # calculate the metrics
    metrics = []
    for i in range(len(all_img)):
        print(str(i) + ' done')
        
        metrics.append(Calc_Metrics(all_img[i], 'all', bm_file, ref_file, 
                                    normal=normal))
        
    # sort the metrics by name of the files:
    names = []
    for i in range(len(all_img)):
        tmp, filename = os.path.split(all_img[i])

        # search for movement type
        test = filename.find('pmcoff')
        if test > 0:
            descr = 'pmcoff_'
        else:
            test = filename.find('pmcon')
            if test > 0:
                descr = 'pmcon_'
            else:
                test = filename.find('MOCO_RETRO_')
                if test > 0:
                    descr = 'MOCO_RETRO_'
                else:
                    print('ERROR: Filename ' + filename + ' does not contain MOCO description')

        # test = filename.find('_RR_')
        test = filename.find('rec-wore')
        if test > 0:
            descr = descr + 'rec-wore_'

        test = filename.find('rec-wre')
        if test > 0:
            descr = descr + 'rec-wre_'

        test = filename.find('run-02')
        if test > 0:
            descr = descr + 'run-02'
        else:
            test = filename.find('run-01')
            if test > 0:
                descr = descr + 'run-01'
            else:
                test = filename.find('run-03')
                if test > 0:
                    descr = descr + 'run-03'
                else:
                    test = filename.find('_SHIFT_')
                    if test > 0:
                        descr = descr + 'SHIFT'                            
                    else:
                        print('EROOR: Filename ' + filename + ' does not contain motion description')

        
        names.append(descr)
    print('names', names)
    names = np.array(names)
    metrics = np.array(metrics)

    ind = np.argsort(names)
    metrics = metrics[ind][::-1]
    names = np.sort(names)[::-1]
    
    a,b = np.shape(metrics)
    save_arr = np.zeros((a,b+1)).astype(str)
    save_arr[:,0]=names
    save_arr[:,1:] = metrics
    header = 'Values of the metrics SSIM, PSNR, Tenengrad, AES, GradEntropy, ImgEntropy, CoEnt for acquisitions specified in first column'
    np.savetxt(out_dir+'Values_Metrics_'+save+tag+'.txt', save_arr, 
               fmt='%s', header=header)
                
    return 0
