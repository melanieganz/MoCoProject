import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
import glob                                                 # Finds all pathnames matching a specified pattern
import nibabel as nib                                       # Package which supports neuroimaging file formats
import shutil                                               # Operations of files
plt.style.use('ggplot')
from recon_register import transformMRI, applyTransformMRI  # Python file with specific functions
                                                            # transformMRI calculates the registration transform
                                                            # applytransformMRI Applies transforms to the scans as well as to the brainmasks
from img_quality_metrics import Comp_All_Metrics            # Python file with specific functions
                                                            # Comp_All_Metrics computes the metrics SSIM, PSNR, Tenengrad, AES, GradEntropy, ImgEntropy, CoEnt for all nifti files specified in all_img.

import fnmatch


def FullAnalysis(sub, nifti_dir, bm_dir, reg_dir, SUBJECTS_DIR, outDir, outDirMetrics, save, 
                 recon_all, register, apply_transform, apply_transform_bm, metrics,show_bm_reg):
    
    '''
    Performs all steps of analysis on prospectively corrected and uncorrected
    scans (if convert, recon_all, transform and metrics are all True)
    
    Parameters
    ----------
    sub : str
        subject ID.
    nfiti_dir : str
        path to nifti directory.
    bm_dir : str
        path to brainmask directory.
    reg_dir : str
        path to directory with registration transforms.
    SUBJECTS_DIR : str
        path to the FreeSurfer subject's directory.
    outDir : str
        path to the output directory for registered scans.
    outDirMetrics : str
        path to the output directory for the metrics.
    save : str
        string appended to the filename of the output, as e.g. date _05_11
    recon_all : bool
        whether to run the still MPRAGE without MoCo though FreeSurfer's 
        recon_all.
    register : bool
        whether to register the images to their respective ground truth.
    apply_transform : bool
        whether to apply the registration transform to the images.
    apply_transform_bm : bool
        whether the registration transform should also be applied to the 
        images. This is not necessary for calculating metrics. Only if one 
        starts from scratch with recon-all etc. for reproducing our calculations.
    metrics : bool
        whether to calculate metrics.
    show_bm_reg : bool
        whether to show and not only save figures that display the registration 
        of the binarized brainmasks with the different sequences
        

    Returns
    -------
    int
        returns 0, when completed.
    
    '''
    
    
    '''Recon_all:'''
    # Run recon all:
    if recon_all:
        list_sequ = os.listdir(nifti_dir) #returns a list containing the names of the entries in the directory given by path
        sequences = [x for x in list_sequ if x.startswith('TCLMOCO')]   #THIS PROBABLY IS NOT UP TO DATE
        sequences = [x for x in list_sequ if x.endswith('.nii')]
        for seq in sequences:
            if seq.startswith('TCLMOCO_OFF_STILL_T1_MPR_3D_SAG'):
                nifti_path = seq 
        
        
        subj_id = sub
        print(subj_id, nifti_path)
        subprocess.run('recon-all -i ' + nifti_path + ' -s ' + subj_id + ' -sd '+SUBJECTS_DIR+' -all -parallel', 
                       shell=True)
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//                  ReconAll done                   //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')
    
    
    
    '''Registration:''' 
    
    if register:
        # get all transforms:
        transformMRI(sub, nifti_dir, outDir)
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//                TransformMRI done                 //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')
        
        
    if apply_transform:
    #    # brainMask all images:
    #    if os.path.exists(bm_dir + '/mri/brainmask_edit.mgz'):
    #        brainmask = bm_dir + '/mri/brainmask_edit.mgz'
    #    else:
    #        brainmask = bm_dir + '/mri/brainmask.mgz'
            
        # if you start from scratch and do not want to use the registration 
        # transforms we provide, the following function needs to be changed
        # to load different registration files!
        print('before')
        applyTransformMRI(nifti_dir, reg_dir, bm_dir, outDir,
                          apply_transform_bm)
        print('after')
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//             applyTransformMRI done               //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')
        
        
    # check brainmasks correctly registered:
    plt.figure(figsize=(10,7))
    i = 1
    #for tag in ['mprage' ]:   #updated name for 1 scan for sub1 to test. should be for all types of scans
    for tag in ['mprage', 'flair', 't2tse', 't1tirm', 't2star']:
        #print('inside mprage bm')
        #print(outDir, tag)
        print(glob.glob(outDir+'*'+tag+'*run-01*'+'*moved.nii'))
        if len(glob.glob(outDir+'*'+tag+'*run-01*'+'*moved.nii')) > 0:
            #print('inside if')
            img_file = glob.glob(outDir+'*'+tag+'*run-01*'+'*moved.nii')[0]
            img = nib.load(img_file).get_fdata().astype(np.uint16)
            bm_file = glob.glob(bm_dir+'*'+tag+'*run-01*'+'mask.nii')[0]
            bm = nib.load(bm_file).get_fdata().astype(np.uint16)


            sh = np.shape(img) #dimensions of image
            sl_dir = np.argmin(sh) #here we figure out which axis/dimension is smallest


            slice_dim = sl_dir
            indx = [slice(None)]
            indx[slice_dim] = int(sh[sl_dir]/2)

            
            plt.subplot(2,4,i) #create plot for each patient
            plt.imshow(np.squeeze(img[indx]), cmap='gray')  #np.squeeze makes us go from shape (1, 256, 256) to shape (256, 256)
            plt.imshow(np.squeeze(bm[indx]), cmap='Reds', alpha=0.4)  #--
            plt.axis('off')
            plt.title('Brainmask '+tag, fontsize = 10)
            i+=1
        
    plt.savefig(outDir+'Check_Brainmask', dpi=300)
    if show_bm_reg:
        plt.show()
    
    
    '''Calculate Metrics:'''   
    if metrics:
        all_img = glob.glob(outDir +'*_moved.nii')
        
        if not os.path.exists(outDirMetrics):
            os.makedirs(outDirMetrics)
        
        # sort different sequences/contrasts:
        test = 0
        if len(all_img) == 0:
            print('No moved images in directory. Check subdirectories.')
        for i in range(len(all_img)):
            tmp, filename = os.path.split(all_img[i])
            # search for different sequence types:
            for tag in ['mprage', 't2tse', 't1tirm', 'flair', 't2star']: #Should be changed to new names
                test = filename.find(tag)
                if test > 0:
                    newpath = tmp + '/' + tag #Create folders if they don't already exist
                    if not os.path.exists(newpath):
                        os.makedirs(newpath)
                    shutil.move(all_img[i], tmp + '/' + tag + '/' + filename) #move images to correct directory
                    test = 1
            if test == 0:
                print('ERROR: Filename ' + filename + ' does not contain the right sequence type description')
        
        
        for tag in ['mprage', 't2tse', 't1tirm', 'flair', 't2star']:
            if os.path.isdir(outDir + tag):
                print('######## ' + tag)
                img_mov = glob.glob(outDir + tag + '/' +'*_moved.nii')
                check = 0
                for file in img_mov:
                    #test = file.find('TCLMOCO_OFF_STILL_')   #Update to new name
                    test = all(fnmatch.fnmatch(file, f"*{keyword}*") for keyword in ['pmcoff', 'run-01'])
                    if test > 0:
                        ref_file = file
                        check = 1
                if check == 0:
                    print('ERROR: no reference Scan in ' + tag)
        
                #bm_mov = glob.glob(bm_dir+'bm_mov_*'+tag+'*.nii')[0]
                bm_mov = glob.glob(bm_dir + '*' + tag + '*.nii')[0]

                Comp_All_Metrics(img_mov, bm_mov, ref_file, outDirMetrics, save, tag)
        
    return 0
       



''' (1) Run analysis on prospectively corrected and uncorrected data:'''
# define specific input parameters for the current run:

#subjs = []
#for i in range(1,10):
#    subjs.append('Subject_0'+str(i))
#for i in range(10,20):
#    subjs.append('Subject_'+str(i))
#for i in range(20,23):
#    subjs.append('Subject_'+str(i))

# DEBUG - for now only run with one subject
sub = 'sub-01'
#Paths should be changed into a loop once we run for all subjects 
root = '/mnt/mocodata1/MoCoHealthy/Public/BIDS/BIDSdata/'
# Path on Windows laptop '//pmod.nru.dk/mocodata1/MoCoHealthy/Public/BIDS/BIDSdata/'
# Path on Unix laptop '/home/melanie/Data/ds004332-download/'

nifti_dir = root+sub+'/anat/'
bm_dir = root+'derivatives/freesurfer/'+sub+'/anat/'
reg_dir = root+'derivatives/freesurfer/'+sub+'/transforms/'
SUBJECTS_DIR = root +'derivatives/freesurfer/'+sub+'/transforms/'
outDir = root+'derivatives/results/registrations/'+sub+'/'
outDirMetrics = root+'derivatives/results/metricsresults/'+sub+'/'

# Original folder structure
#for sub in subjs:
    #nifti_dir = '../BIDSdata_defaced/'+sub+'/'
    #bm_dir = '../Brainmasks/'+sub+'/'
    #reg_dir = '../RegistrationTransforms/'+sub+'/'
    #SUBJECTS_DIR = '../Results/Data_Recon_All/'
    #outDir = '../Results/Registrations/'+sub+'/'
    #outDirMetrics = '../Results/Metrics_Results/'+sub+'/'

save = '2023_06_08'
    
    # Which steps to perform:
    # by default we have performed registrations in FreeSurfer and provide those for you
    # if wanted one could set recon_all and register to True and redo this analysis as well,
    # but this would then entail running recon-all on all Still images and hence take a while

    #Variables before subject

recon_all = False
register = False
apply_transform_bm = False
    
# steps to basically always perform
apply_transform = True
metrics = True
show_bm_reg = False

#Run the function
FullAnalysis(sub, nifti_dir, bm_dir, reg_dir, SUBJECTS_DIR, outDir, outDirMetrics, save, recon_all=recon_all, register=register, apply_transform=apply_transform, apply_transform_bm=apply_transform_bm, metrics=metrics, show_bm_reg=show_bm_reg)



