import os
import glob
import subprocess
import numpy as np
import nibabel as nib


def transformMRI(sub, niftiDir, outDir):
    '''
    Calculates the registration transform: using FreeSurfer functions 
    bbregister for MPRAGE and still MOCO_OFF of remaining scans and 
    robust register for remaining scans.

    Parameters
    ----------
    sub : str
        subject ID.
    niftiDir : str
        path to Nifti directory.
    outDir : str
        path to output directory.

    Returns
    -------
    int
        returns 0, when completed.

    '''

    # use bbregister for T1 MPRAGE scans and MOCO_OFF_STILL_ scans 
    # corresponding to remaining sequences (in order to later apply to 
    # brainmask):
    regDir = outDir + 'regs/'
    if not os.path.exists(regDir):
        print('New regDir created')
        os.makedirs(regDir)
    files = glob.glob(niftiDir+"*T1_MPR_*.nii")

    for i in range(0, len(files)):
        movImg = niftiDir + os.path.basename(files[i])
        regname = regDir + os.path.basename(files[i]) + '.lta'

        subprocess.run('bbregister --s ' + sub + ' --mov '+  movImg + ' --reg ' + regname + ' --t1 --init-best-header', shell=True)
        print(i, ' done')

    print('##################')
    print('bbregister for T1_MPR_ done')
    print('##################')

    for tag in ['T2_TSE_', 'T1_TIRM_', 'T2_FLAIR', 'T2STAR', 'EPI_SWI', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
        if len(glob.glob(niftiDir+"*MOCO_OFF_STILL_*"+tag+"*.nii"))>0:
            still_img = glob.glob(niftiDir+"*MOCO_OFF_STILL_*"+tag+"*.nii")[0]
            movImg = niftiDir + os.path.basename(still_img)
            regname = regDir + os.path.basename(still_img) + '.lta'

            if tag in ['T1_TIRM_']:
                subprocess.run('bbregister --s ' + sub + ' --mov '+  movImg + ' --reg ' + regname + ' --t1 --init-best-header', shell=True)
            else:
                subprocess.run('bbregister --s ' + sub + ' --mov '+  movImg + ' --reg ' + regname + ' --t2 --init-best-header', shell=True)

    print('##################')
    print('bbregister for MOCO_OFF_STILL_ of remaining scans done')
    print('##################')


    # use robust register for the remaining sequences:
    regDir = outDir + 'regs_robust/'
    if not os.path.exists(regDir):
        print('New regDir created')
        os.makedirs(regDir)

    for tag in ['T2_TSE_', 'T1_TIRM_', 'T2_FLAIR', 'T2STAR', 'EPI_SWI', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
        if len(glob.glob(niftiDir+"*"+tag+"*.nii"))>0:
            files = glob.glob(niftiDir+"*"+tag+"*.nii")
            targImg = glob.glob(niftiDir+"*MOCO_OFF_STILL_*"+tag+"*.nii")[0]

            for i in range(0, len(files)):
                movImg = niftiDir + os.path.basename(files[i])
                regname = regDir + os.path.basename(files[i]) + '.lta'

                subprocess.run('mri_robust_register --mov '+  movImg + ' --dst ' + targImg + ' --lta ' + regname + ' --satit', shell=True)
                print(i, ' done')

            print('##################')
            print('robust register for ', tag, ' done')
            print('##################')

    return 0


def applyTransformMRI(niftiDir, brainmask, outDir):
    '''
    Applies transforms (saved in outDir/regs or outDir/regs_robust) to the 
    scans as well as to the brainmasks

    Parameters
    ----------
    niftiDir : str
        directory where nifti images are stored.
    brainmask : str
        directory where brainmask for reference image is stored.
    outDir : str
        directory where moved and masked images will be saved.

    Returns
    -------
    int
        returns 0, when completed.

    '''
    
    # binarize brainmask and transform T1_MPR:
    # output: brainmask_bin.nii and *T1_MPR_*_moved.nii
    regDir = outDir + 'regs/'
    subDir, tail = os.path.split(brainmask)
    name, ext = os.path.splitext(tail)
    brainmask_bin_nii = outDir + 'brainmask_bin.nii'
    T1 = subDir + '/T1.mgz'

    # binarize brainmask:
    subprocess.run('mri_binarize --i ' + brainmask + ' --o ' + brainmask_bin_nii + ' --match 0 --inv', shell=True)

    files = glob.glob(niftiDir+"*T1_MPR_*.nii")

    for i in range(0, len(files)):
        vol = niftiDir+os.path.basename(files[i])
        name2, ext2 = os.path.splitext(os.path.basename(files[i]))
        vol_moved = outDir + name2 + '_moved' + ext2
        regname = regDir + os.path.basename(files[i]) + '.lta'

        # transform:
        subprocess.run('mri_vol2vol --mov ' + vol + ' --targ ' + T1 + ' --o ' + vol_moved + ' --lta ' + regname, shell=True)

        print(i, ' done')

    print('##################')
    print('Transform applied for T1_MPR_')
    print('##################')


    # transform brainmask into still/MoCo_off domain of remaining scans:
    for tag in ['T2_TSE_', 'T1_TIRM_', 'T2_FLAIR', 'T2STAR', 'EPI_SWI', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
        if len(glob.glob(niftiDir+"*MOCO_OFF_STILL_*"+tag+"*.nii"))>0:
            still_img = glob.glob(niftiDir+"*MOCO_OFF_STILL_*"+tag+"*.nii")[0]
            bm_mov = outDir + 'bm_mov_'+os.path.basename(still_img)
            T2 = niftiDir + os.path.basename(still_img)
            regname = regDir + os.path.basename(still_img) + '.lta'
            name2, ext2 = os.path.splitext(os.path.basename(still_img))

            # transform:
            subprocess.run('mri_vol2vol --mov ' + T2 + ' --targ ' + brainmask_bin_nii + ' --o ' + bm_mov + ' --lta ' + regname + ' --inv --nearest', shell=True)

            print('##################')
            print('Brainmask transformed for ', tag)
            print('##################')


    # transform remaining scans:
    # output: *tag*_moved.nii
    regDir = outDir + 'regs_robust/'

    for tag in ['T2_TSE_', 'T1_TIRM_', 'T2_FLAIR', 'T2STAR', 'EPI_SWI', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
        if len(glob.glob(niftiDir+"*MOCO_OFF_STILL_*"+tag+"*.nii"))>0:
            targImg = glob.glob(niftiDir+"*MOCO_OFF_STILL_*"+tag+"*.nii")[0]
            files = glob.glob(niftiDir+"*"+tag+"*.nii")

            for i in range(0, len(files)):
                vol = niftiDir+os.path.basename(files[i])
                name2, ext2 = os.path.splitext(os.path.basename(files[i]))
                vol_moved = outDir + name2 + '_moved' + ext2
                regname = regDir + os.path.basename(files[i]) + '.lta'

                # transform:
                subprocess.run('mri_vol2vol --mov ' + vol + ' --targ ' + targImg + ' --o ' + vol_moved + ' --lta ' + regname, shell=True)

                print(i, ' done')

    return 0



def transformMRIRetro(sub, niftiDir, outDir, corr=True):
    '''
    Calculates the registration transform for retrospectively corrected data:
    using FreeSurfer functions bbregister for MPRAGE and still MOCO_OFF scans 
    of remaining sequences and robust register for remaining scans.

    Parameters
    ----------
    sub : str
        subject ID.
    niftiDir : str
        path to Nifti directory.
    outDir : str
        path to output directory.
    corr : bool (optional)
        if it is True, then regnames are adjusted to deal with the fact that
        retrospectively reconstructed scans have a different filename after
        bias field correction.

    Returns
    -------
    int
        returns 0, when completed.

    '''

    # bbregister for T1 MPRAGE and MOCO_OFF_STILL_ FLAIR scan:
    regDir = outDir + 'regs_retro/'
    if not os.path.exists(regDir):
        print('New regDir created')
        os.makedirs(regDir)
    files = glob.glob(niftiDir+"*T1_MPR_*.nii")

    for i in range(0, len(files)):
        movImg = niftiDir + os.path.basename(files[i])
        regname = regDir + os.path.basename(files[i]) + '.lta'
        if corr == True:
            regname = regDir + os.path.basename(files[i])[14:] + '.lta'

        subprocess.run('bbregister --s ' + sub + ' --mov '+  movImg + ' --reg ' + regname + ' --t1 --init-best-header', shell=True)
        print(i, ' done')

    print('##################')
    print('bbregister for T1_MPR_Retro done')
    print('##################')


    # robust register for all T2 FLAIR scans:
    regDir = outDir + 'regs_robust_retro/'
    if not os.path.exists(regDir):
        print('New regDir created')
        os.makedirs(regDir)


    if len(glob.glob(niftiDir+"*T2_FLAIR*.nii"))>0:
        files = glob.glob(niftiDir+"*T2_FLAIR*.nii")
        targImg = glob.glob(niftiDir+"*MOCO_OFF_STILL_*T2_FLAIR*.nii")[0]

        for i in range(0, len(files)):
            movImg = niftiDir + os.path.basename(files[i])
            regname = regDir + os.path.basename(files[i]) + '.lta'
            if corr == True:
                regname = regDir + os.path.basename(files[i])[14:] + '.lta'

            subprocess.run('mri_robust_register --mov '+  movImg + ' --dst ' + targImg + ' --lta ' + regname + ' --iscale --satit', shell=True)
            print(i, ' done')

        print('##################')
        print('robust register for T2_FLAIR_Retro done')
        print('##################')

    return 0



def applyTransformMRIRetro(niftiDir, brainmask, outDir, corr=True):
    '''
    Applies transforms for retrospectively corrected scans (saved in
    outDir/regs or outDir/regs_robust) and transforms the brainmasks

    Parameters
    ----------
    niftiDir : str
        directory where nifti images are stored.
    brainmask : str
        directory where brainmask for reference image is stored.
    outDir : str
        directory where moved and masked images will be saved.
    corr : bool (optional)
        if it is True, then regnames are adjusted to deal with the fact that
        retrospectively reconstructed scans have a different filename after
        bias field correction.

    Returns
    -------
    int
        returns 0, when completed.

    '''
    # binarize brainmask and transform T1_MPR
    # --> output: brainmask_bin.nii and *T1_MPR_*_moved.nii
    regDir = outDir + 'regs_retro/'
    subDir, tail = os.path.split(brainmask)
    name, ext = os.path.splitext(tail)
    brainmask_bin_nii = outDir + name + '_bin.nii'
    T1 = subDir + '/T1.mgz'

    # binarize brainmask:
    subprocess.run('mri_binarize --i ' + brainmask + ' --o ' + brainmask_bin_nii + ' --match 0 --inv', shell=True)

    files = glob.glob(niftiDir+"*T1_MPR_*.nii")

    for i in range(0, len(files)):
        vol = niftiDir+os.path.basename(files[i])
        name2, ext2 = os.path.splitext(os.path.basename(files[i]))
        vol_moved = outDir + 'retro_scans/' + name2 + '_moved' + ext2
        regname = regDir + os.path.basename(files[i]) + '.lta'
        if corr == True:
            regname = regDir + os.path.basename(files[i])[14:] + '.lta'

        # transform:
        subprocess.run('mri_vol2vol --mov ' + vol + ' --targ ' + T1 + ' --o ' + vol_moved + ' --lta ' + regname, shell=True)

        print(i, ' done')

    print('##################')
    print('Transform applied for T1_MPR_')
    print('##################')



    # now move T2 FLAIR scans:
    regDir = outDir + 'regs_robust_retro/'

    if len(glob.glob(niftiDir+"*MOCO_OFF_STILL_*T2_FLAIR*.nii"))>0:
        targImg = glob.glob(niftiDir+"*MOCO_OFF_STILL_*T2_FLAIR*.nii")[0]
        files = glob.glob(niftiDir+"*T2_FLAIR*.nii")

        for i in range(0, len(files)):
            vol = niftiDir+os.path.basename(files[i])
            name2, ext2 = os.path.splitext(os.path.basename(files[i]))
            vol_moved = outDir + 'retro_scans/' + name2 + '_moved' + ext2
            regname = regDir + os.path.basename(files[i]) + '.lta'
            if corr == True:
                regname = regDir + os.path.basename(files[i])[14:] + '.lta' #check that working

            # transform:
            subprocess.run('mri_vol2vol --mov ' + vol + ' --targ ' + targImg + ' --o ' + vol_moved + ' --lta ' + regname, shell=True)

            print(i, ' done')
    
        # use the brainmask transform calculated on data reconstructed by the 
        # scanner, but flip it in all directin so that it matches 
        # retrospectively reconstructed FLAIR scans):
        if corr == False:
            bm_file = glob.glob(outDir+'bm_mov_*TCLMOCO_OFF_STILL_T2_FLAIR*.nii')[0]
            bm = nib.load(bm_file).get_fdata().astype(np.uint16)
            bm_tr = bm[::-1,::-1,::-1]
            np.save(outDir+'bm_mov_retro_TCLMOCO_OFF_STILL_T2_FLAIR', bm_tr)
        
    print('##############')
    print('Transform applied for T2 FLAIR')
    print('##############')

    return 0



def Run_Recon_All_Again(fails):
    for f in fails:
        subprocess.run('recon-all -s ' + f + ' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all  -parallel', shell=True)   # -sd /data1/hannah/Old_Setup/Data_Recon_All/

    return 0



def Run_Long_Stream(name):

    # Run base recon all to create whithin-subject template:
    #  MoCo OFF scans
    tp = '-tp '+name # MoCo OFF Still has different name
    for mov in ['NOD_RR_', 'SHAKE_RR_']:
        tp += ' -tp '+'X_OFF_'+ mov+name
    for mov in ['STILL_', 'NOD_RR_', 'SHAKE_RR_']:
        tp += ' -tp '+'X_ON_'+ mov+name

    if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/X_ON_STILL_Vol_10/Data_Recon_All/X_BASE_OFF_'+name)==False:
        subprocess.run('recon-all -base X_BASE_'+name+' '+tp+' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all  -parallel', shell=True)


    # Longitudinal runs:
    if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/'+name+'.long.X_BASE_OFF_'+name)==False:
        subprocess.run('recon-all -long '+name+' X_BASE_'+name+' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all  -parallel', shell=True) # MoCo OFF Still has different name
    for mov in ['NOD_RR_', 'SHAKE_RR_']:
        if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/X_OFF_'+mov+name+'.long.X_BASE_OFF_'+name)==False:
            subprocess.run('recon-all -long X_OFF_'+mov+name+' X_BASE_'+name+' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all  -parallel', shell=True)
    for mov in ['STILL_', 'NOD_RR_', 'SHAKE_RR_']:
        if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ON_'+mov+name+'.long.X_BASE_ON_'+name)==False:
            subprocess.run('recon-all -long X_ON_'+mov+name+' X_BASE_'+name+' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all  -parallel', shell=True)

    return 0














