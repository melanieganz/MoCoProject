'''
Script for evaluating image quality metrics on the challenge test submissions
'''

import numpy as np
import subprocess
import os
from Utils import RegisterNifti, Comp_Metrics


# define directories etc:
subm_dir = 'Submissions'    # path to the folder containing subfolders for each
                            # participant.
intern_dir = 'Internal'     # path to the folder containing ground truth (still)
                            # scans for each participant
out_dir = 'Processed_Submissions'   # path to the folder where all intermediate
                                    # files for processing should be saved
participants = ['Team1', 'Team2']   # list of team names
test_subj = []      # list of IDs for test subjects
for i in range(1,10):
    test_subj.append('Test_0'+str(i))
for i in range(10,21):
    test_subj.append('Test_'+str(i))
    

# specify whether the reconstruction should be registered to the ground
# truth scan - only needs to be done for test data set:
register_niftis = False
# specify whether the registration should be manually checked - for that
# the moved image and ground truth scan will be loaded with Freeview
check_registr = False


for part in participants:
    # check if submissions complete:
    missing_t1, missing_t2 = [], []
    print('Missing reconstructions for '+part+':')
    for sub in test_subj:
        for task, missing in zip(['t1', 't2'], [missing_t1, missing_t2]):
            if not os.path.exists(subm_dir+part+'/'+sub+'/Recon_nod_'+task+'.nii'):
                missing.append(sub+'/Recon_nod_'+task+'.nii')
                print(sub+'/Recon_nod_'+task+'.nii')

            else:
                # register the submission to the ground truth scan (optional)
                if register_niftis:
                    RegisterNifti(sub, task, subm_dir+part+'/', intern_dir,
                                  out_dir+part+'/')

    # manually check the registrations:
    if check_registr:
        for sub in test_subj:
            for task in ['t1', 't2']:
                if sub+'/Recon_nod_'+task+'.nii' not in missing_t1+missing_t2:
                    mov = out_dir+part+'/'+sub+'/Recon_nod_'+task+'_mov.nii'
                    ref = intern_dir+sub+'/Scan_still_'+task+'_off.nii'
                    subprocess.run('freeview -v '+mov+' -v '+ref, shell=True)


    # calculate SSIM for ranking, as well as PSNR and Tenengrad for additional
    # analysis (not considered in the evaluation/ranking process):
    for task in ['t1', 't2']:
        ssim, psnr, tg = [], [], []
        for sub in test_subj:
            if sub+'/Recon_nod_'+task+'.nii' not in missing_t1+missing_t2:
                img = out_dir+part+'/'+sub+'/Recon_nod_'+task+'_mov.nii'
                ref = intern_dir+sub+'/Scan_still_'+task+'_off.nii'
                bm = intern_dir+sub+'/Brainmask_still_'+task+'_off.nii'
                s, p, t = Comp_Metrics(img, ref, bm)
                ssim.append(s)
                psnr.append(p)
                tg.append(t)
            else:
                # if no submission for that subject, this is punished with a
                # value of 0 for that subject:
                ssim.append(0)
                psnr.append(0)
                tg.append(0)

                
        # calculate average metric values and store them in footer of the text
        # file:
        footer = 'Average over all subjects: '
        for m in [ssim, psnr, tg]:
            footer += str(np.mean(m))+' '
        
        
        # save the calculated metrics in a text file:
        save_arr = np.array([test_subj, ssim, psnr, tg]).T
        np.savetxt(out_dir+part+'/Values_'+task+'.txt', save_arr, fmt='%s', 
                   header='Team '+part+'\n Values of the metrics SSIM, PSNR and TG for subjects specified in first column',
                   footer=footer)
