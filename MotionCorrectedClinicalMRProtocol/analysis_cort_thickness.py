'''Warning - this script reruns FreeSurfers recon-all on all images and therefore takes very long to run'''
import numpy as np
import subprocess
import os
import datetime
import glob
from recon_register import Run_Recon_All_Again, Run_Long_Stream
from utils import SortFiles

run_RR = False
run_no_RR = False
run_long = False

names = []
for i in range(1,10):
    names.append('Subject_0'+str(i)+'/')
for i in range(10,20):
    names.append('Subject_'+str(i)+'/')
for i in range(20,23):
    names.append('Subject_'+str(i)+'/')
    
    
''' (1)  Run ReconAll for the motion RR scans in order to compare freesurfer estimates'''

if run_RR == True:
    for name in names:
        # define directories: 
        nifti_dir = '../BIDSdata_defaced/'+name+'/'
        
        
            
        list_sequ = os.listdir(nifti_dir)
        sequences = [x for x in list_sequ if x.startswith('TCLMOCO')] 
        
        # check that nifti directory exists, otherwise make a new nifti directory 
        # for this subject: 
        if not os.path.exists(nifti_dir):
            os.makedirs(nifti_dir)
        
        
        for seq in sequences:
            if 'T1_MPR' not in seq:
                continue
            #dcm = os.listdir(nifti_dir + seq)[0]
            #in_volume = dicom_dir + seq + '/' + dcm
            out_volume = nifti_dir + 'TCL'+name+'_' + seq + '.nii'
            if seq.startswith('TCLMOCO_OFF_STILL_T1_MPR'):
                continue
            
        
            # Run recon all:
            descr = ['ON_STILL_', 'ON_NOD_', 'ON_SHAKE_', 'OFF_NOD_', 'OFF_SHAKE_']
            for d in descr:        
                if d+'T1_MPR' in seq:
                    if d == 'ON_STILL_':
                        subj_id = 'X_'+d+name
                    else:
                        if 'RR' in seq:
                            subj_id = 'X_'+d+'RR_'+name
                        else:
                            print('No RR scan for ', d)
                            continue
        
                    print(subj_id)
                    if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/'+subj_id):
                        continue
                
                    subprocess.run('recon-all -i ' + out_volume + ' -s ' + subj_id + ' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all -parallel', 
                                   shell=True)
                    
                    print('############')
                    print(name, ' done')
                    print('############')



''' (2) Run ReconAll for the motion MoCo ON/OFF scans with REAC in order to compare freesurfer estimates'''

if run_no_RR:
    for name in names:
        # define directories:
        #dicom_dir = '/mnt/mocodata1/Data_Analysis/DICOMS/'+name+'/'   
        nifti_dir = '../BIDSdata_defaced/'+name+'/'
        
        
            
        list_sequ = os.listdir(nifti_dir)
        sequences = [x for x in list_sequ if x.startswith('TCLMOCO')] 
        
        # check that nifti directory exists, otherwise make a new nifti directory 
        # for this subject: 
        if not os.path.exists(nifti_dir):
            os.makedirs(nifti_dir)
        
        
        for seq in sequences:
            if 'T1_MPR' not in seq:
                continue
            #dcm = os.listdir(dicom_dir + seq)[0]
            #in_volume = dicom_dir + seq + '/' + dcm
            out_volume = nifti_dir +'TCL'+name[:-1]+'_' + seq + '.nii'
            if seq.startswith('TCLMOCO_OFF_STILL_T1_MPR'):
                continue
            
        
            # Run recon all:
            descr = ['OFF_NOD_', 'OFF_SHAKE_']
            
            for d in descr:
                if d in seq and 'RR' not in seq:
                    subj_id = 'X_'+d+name
                    
                    print(subj_id, out_volume)
            
                    if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/'+subj_id):
                        continue
                
                    subprocess.run('recon-all -i ' + out_volume + ' -s ' + subj_id + ' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all -parallel', 
                                   shell=True)
                    
                    print('############')
                    print(name, ' done')
                    print('############')
                    
                    
                    with open('/mnt/mocodata1/Data_Analysis/Surface_Estimates/Status_Recon_All.txt', 'a') as f:
                        f.write(subj_id+' done at '+str(datetime.datetime.now())+'\n')



''' (3) Check whether cross-sectional runs were successful. If not, run them 
again. Otherwise run longitudinal stream. Afterwards, check longitudinal runs '''

outDir = '/mnt/mocodata1/Data_Analysis/Surface_Estimates/'

if run_long:
    for name in names:
        descr = ['ON_STILL_', 'ON_NOD_',  'OFF_NOD_', 'ON_SHAKE_', 'OFF_SHAKE_']
    
        # Test whether cross recon - all successful:
        fails = []
    
        for d in descr:
            if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/X_'+d+name+'/scripts/recon-all.error'):
                print('Recon all failed for ', 'X_', d, name)
                fails.append('X_'+d+name)
            if d != 'ON_STILL_':
                if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/X_'+d+'RR_'+name+'/scripts/recon-all.error'):
                    print('Recon all failed for ', 'X_'+d+'RR_'+name)
                    fails.append('X_'+d+'RR_'+name)
    
        # run recon All again for failed runs:
        Run_Recon_All_Again(fails)
    
    
        # run longitudinal processing stream:
        Run_Long_Stream(name)
    
        # inspect outcomes of base and long!!!!!


    # Check longitudinal runs:
    for name in names:
        descr = ['ON_STILL_', 'ON_NOD_',  'OFF_NOD_', 'ON_SHAKE_', 'OFF_SHAKE_']
    
        fails = []
        for d in descr:
            if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/X_'+d+name[:-1]+'.long.X_BASE_'+name[:-1]+'/scripts/recon-all.error'):
                print('Recon all failed for ', 'X_', d, name, '.long')
                fails.append('X_'+d+name[:-1])
            if d != 'ON_STILL_':
                if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/X_'+d+'RR_'+name[:-1]+'.long.X_BASE_'+name[:-1]+'/scripts/recon-all.error'):
                    print('Recon all failed for ', 'X_'+d+'RR_'+name+'.long')
                    fails.append('X_'+d+'RR_'+name[:-1])
    
        if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/'+name[:-1]+'.long.X_BASE_'+name[:-1]+'/scripts/recon-all.error'):
            print('Recon all failed for '+name[:-1]+'.long.')
            fails.append(name[:-1])
    
        if os.path.exists('/mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/X_BASE_'+name+'scripts/recon-all.error'):
            print('Recon all failed for X_BASE_'+name)
            fails.append('X_BASE_'+name[:-1])
    
        
        # run longitudinal recon-all again without parallel option for failed runs:
        for f in fails:
            if f != 'X_BASE_'+name[:-1]:
                subprocess.run('recon-all -long '+f+' X_BASE_'+name[:-1]+' -sd /mnt/mocodata1/Data_Analysis/Data_Recon_All/Longitudinal/ -all', shell=True)
    
                with open('/mnt/mocodata1/Data_Analysis/Surface_Estimates/Status_Long_ReRun.txt', 'a') as file:
                    file.write(f+'_RR re-run done at '+str(datetime.datetime.now())+'\n')
        
    


''' (4) Calculate average cortical thickness across all volunteers and all 
cortical regions - first for cross-sectional, then for longitudinal runs:'''

# compute only for ground truth scans
for name in names:
    subprocess.run('mris_anatomical_stats -log '+outDir+'AnatomicalStats/Stats_'+name[:-1]+'.txt '+name[:-1]+ ' lh', shell=True)
    subprocess.run('mris_anatomical_stats -log '+outDir+'AnatomicalStats/Stats_rh_'+name[:-1]+'.txt '+name[:-1]+ ' rh', shell=True)
    a=1

#now load all data:
all_th = []
for name in names:
    tmp = np.loadtxt(outDir+'AnatomicalStats/Stats_'+name[:-1]+'.txt', skiprows=1, usecols=2, max_rows=1)
    all_th.append(tmp)
    tmp = np.loadtxt(outDir+'AnatomicalStats/Stats_rh_'+name[:-1]+'.txt', skiprows=1, usecols=2, max_rows=1)
    all_th.append(tmp)

print('CROSS RUNS:')
print('Mean thickness: ', np.mean(all_th), ' +- std ', np.std(all_th))
print('20 percent: ', 0.2*np.mean(all_th))
print('4 percent: ', 0.04*np.mean(all_th))


# compute for longitudinal scans:
for name in names:
    subprocess.run('mris_anatomical_stats -log '+outDir+'AnatomicalStats/Stats_'+name[:-1]+'long.txt '+name[:-1]+'.long.X_BASE_'+name[:-1]+ ' lh', shell=True)
    subprocess.run('mris_anatomical_stats -log '+outDir+'AnatomicalStats/Stats_rh_'+name[:-1]+'long.txt '+name[:-1]+'.long.X_BASE_'+name[:-1]+  ' rh', shell=True)
    a=1

#now load all data:
all_th = []
for name in names:
    tmp = np.loadtxt(outDir+'AnatomicalStats/Stats_'+name[:-1]+'long.txt', skiprows=1, usecols=2)
    all_th.append(tmp)
    tmp = np.loadtxt(outDir+'AnatomicalStats/Stats_rh_'+name[:-1]+'long.txt', skiprows=1, usecols=2)
    all_th.append(tmp)

print('LONGITUDINAL RUNS:')
print('Mean thickness: ', np.mean(all_th), ' +- std ', np.std(all_th))
print('20 percent: ', 0.2*np.mean(all_th))
print('4 percent: ', 0.04*np.mean(all_th))




''' (5) GLM analysis cross-sectional runs'''

run = True
if run == True:

    '''Create FSGD file / design matrix for scans without reacquisition: '''
    #RMS for MoCo Off still:
    descr = ['ON_STILL_', 'ON_NOD_',  'OFF_NOD_', 'ON_SHAKE_', 'OFF_SHAKE_']
    motion = ['STILL', 'NOD', 'NOD', 'SHAKE', 'SHAKE']
    files = glob.glob('/mnt/mocodata1/Data_Analysis/Motion_Estimates/Comparison/MotionMetrics_STILL/T1_MPR*.txt')
    file = SortFiles(files)[0]
    metrics_ref = np.loadtxt(file, unpack=True, usecols=0)
    subj = []
    for i in range(1,10):
        subj.append('HC_0'+str(i)+'/')
    for i in range(10,20):
        subj.append('HC_'+str(i)+'/')
    for i in range(20,23):
        subj.append('HC_'+str(i)+'/')
    subj = np.array(subj)

    for d, m in zip(descr, motion):
        lines = []
        files = glob.glob('/mnt/mocodata1/Data_Analysis/Motion_Estimates/Comparison/MotionMetrics_'+m+'/T1_MPR*.txt')
        files = [f for f in files if 'mid' not in f]    # sort out all files 
        # where displacement was caluclated relative to mid of acquisition
        file = SortFiles(files)[0]
        if d[:2] == 'OF':
            metrics = np.loadtxt(file, unpack=True, usecols=0)
        if d[:2] == 'ON':
            metrics = np.loadtxt(file, unpack=True, usecols=3)

        for name in names:
            s = np.where(subj==name)[0]
            RMS_ref = metrics_ref[s]
            RMS = metrics[s]

            lines.append('Input '+name+' Main '+str(RMS_ref[0]))
            if 'STILL' in d:
                lines.append('Input X_'+d+name+' Main '+str(RMS[0]))
            else:
                lines.append('Input X_'+d+'RR_'+name+' Main '+str(RMS[0]))
        print(d + ' done')

        save = np.array(lines)
        head = 'GroupDescriptorFile 1 \nMeasurementName thickness \nClass Main'
        head += ' \nVariables RMS'
        np.savetxt(outDir+'MOCO_'+d+'.fsgd', save, header=head, comments='', fmt='%s')


    '''Create FSGD file / design matrix for scans with reacquisition: '''
    descr = ['ON_NOD_', 'ON_SHAKE_', 'OFF_NOD_', 'OFF_SHAKE_']
    motion = ['NOD', 'SHAKE', 'NOD', 'SHAKE']
    files = glob.glob('/mnt/mocodata1/Data_Analysis/Motion_Estimates/Comparison/MotionMetrics_STILL/T1_MPR*.txt')
    file = SortFiles(files)[0]
    metrics_ref = np.loadtxt(file, unpack=True, usecols=0)
    subj = []
    for i in range(1,10):
        subj.append('HC_0'+str(i)+'/')
    for i in range(10,20):
        subj.append('HC_'+str(i)+'/')
    for i in range(20,23):
        subj.append('HC_'+str(i)+'/')
    subj = np.array(subj)

    for d, m in zip(descr, motion):
        lines = []
        files = glob.glob('/mnt/mocodata1/Data_Analysis/Motion_Estimates/Comparison/MotionMetrics_'+m+'/T1_MPR*.txt')
        file = SortFiles(files)[0]
        if d[:2] == 'OF':
            metrics = np.loadtxt(file, unpack=True, usecols=0)
        if d[:2] == 'ON':
            metrics = np.loadtxt(file, unpack=True, usecols=3)

        for name in names:
            s = np.where(subj==name)[0]
            RMS_ref = metrics_ref[s]
            RMS = metrics[s]

            lines.append('Input '+name+' Main '+str(RMS_ref[0]))
            lines.append('Input X_'+d+name+' Main '+str(RMS[0]))
        print(d + ' done')

        save = np.array(lines)
        head = 'GroupDescriptorFile 1 \nMeasurementName thickness \nClass Main'
        head += ' \nVariables RMS'
        np.savetxt(outDir+'MOCO_'+d+'_REAC_ON.fsgd', save, header=head, comments='', fmt='%s')



    '''Assemble the data: '''
    # 10 to 15 minutes
    descr = ['MOCO_ON_STILL_', 'MOCO_ON_NOD_',  'MOCO_OFF_NOD_', 'MOCO_ON_SHAKE_', 'MOCO_OFF_SHAKE_', 'MOCO_ON_NOD__REAC_ON', 'MOCO_OFF_NOD__REAC_ON', 'MOCO_ON_SHAKE__REAC_ON',  'MOCO_OFF_SHAKE__REAC_ON']
    for d in descr:
        # lh
        # resample to fsaverage and put all in one file:
        subprocess.run('mris_preproc --fsgd '+outDir+d+'.fsgd --target fsaverage --hemi lh --meas thickness --out '+outDir+'lh.'+d+'.thickness.00.mgh', shell=True)
        # smooth each subject's resampled data by 10mm FWHM
        subprocess.run('mri_surf2surf --hemi lh --s fsaverage --sval '+outDir+'lh.'+d+'.thickness.00.mgh --fwhm 10 --cortex --tval '+outDir+'lh.'+d+'.thickness.10.mgh', shell=True)

        # rh
        # resample to fsaverage and put all in one file:
        subprocess.run('mris_preproc --fsgd '+outDir+d+'.fsgd --target fsaverage --hemi rh --meas thickness --out '+outDir+'rh.'+d+'.thickness.00.mgh', shell=True)
        # smooth each subject's resampled data by 10mm FWHM
        subprocess.run('mri_surf2surf --hemi rh --s fsaverage --sval '+outDir+'rh.'+d+'.thickness.00.mgh --fwhm 10 --cortex --tval '+outDir+'rh.'+d+'.thickness.10.mgh', shell=True)



    '''GLM Analysis: '''
    for d in descr:
        print(d)
        # fit glm lh
        subprocess.run('mri_glmfit --y '+outDir+'lh.'+d+'.thickness.10.mgh --fsgd '+outDir+d+'.fsgd --C '+outDir+'lh-Avg-thickness-RMS-Cor.mtx --surf fsaverage lh --cortex --glmdir '+outDir+'lh.'+d+'.glmdir', shell=True)

        # fit glm rh
        subprocess.run('mri_glmfit --y '+outDir+'rh.'+d+'.thickness.10.mgh --fsgd '+outDir+d+'.fsgd --C '+outDir+'rh-Avg-thickness-RMS-Cor.mtx --surf fsaverage rh --cortex --glmdir '+outDir+'rh.'+d+'.glmdir', shell=True)

        # fdr correction (lh and rh together):
        subprocess.run('mri_fdr --i '+outDir+'lh.'+d+'.glmdir/lh-Avg-thickness-RMS-Cor/sig.mgh nomask '+outDir+'lh.'+d+'.glmdir/lh-Avg-thickness-RMS-Cor/sig_fdr.mgh --i '+outDir+'rh.'+d+'.glmdir/rh-Avg-thickness-RMS-Cor/sig.mgh nomask '+outDir+'rh.'+d+'.glmdir/rh-Avg-thickness-RMS-Cor/sig_fdr.mgh --fdr 0.05', shell=True)

        # view the results (significance maps and beta - change)- opens new freeview window for each!!!
        folder = outDir+'lh.'+d+'.glmdir/lh-Avg-thickness-RMS-Cor/'
        subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'sig_fdr.mgh:overlay_threshold=1.301,6 -viewport 3d', shell=True)
        folder = outDir+'lh.'+d+'.glmdir/'
        subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'beta.mgh:overlay_threshold=0.04,0.2 -viewport 3d', shell=True)

        # rh
        folder = outDir+'rh.'+d+'.glmdir/rh-Avg-thickness-RMS-Cor/'
        #subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/rh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'sig_fdr.mgh:overlay_threshold=1.301,6 -viewport 3d', shell=True)
        folder = outDir+'rh.'+d+'.glmdir/'
        #subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/rh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'beta.mgh:overlay_threshold=0.04,0.2 -viewport 3d', shell=True)



'''  (6) GLM analysis longitudinal runs '''

run = False
if run == True:

    '''Create FSGD file / design matrix for scans without reacquisition: '''
    #RMS for MoCo Off still:
    descr = ['ON_STILL_', 'ON_NOD_',  'OFF_NOD_', 'ON_SHAKE_', 'OFF_SHAKE_']
    motion = ['STILL', 'NOD', 'NOD', 'SHAKE', 'SHAKE']
    files = glob.glob('/mnt/mocodata1/Data_Analysis/Motion_Estimates/Comparison/MotionMetrics_STILL/T1_MPR*.txt')
    file = SortFiles(files)[0]
    metrics_ref = np.loadtxt(file, unpack=True, usecols=0)
    subj = []
    for i in range(1,10):
        subj.append('HC_0'+str(i)+'/')
    for i in range(10,20):
        subj.append('HC_'+str(i)+'/')
    for i in range(20,23):
        subj.append('HC_'+str(i)+'/')
    subj = np.array(subj)

    for d, m in zip(descr, motion):
        lines = []
        files = glob.glob('/mnt/mocodata1/Data_Analysis/Motion_Estimates/Comparison/MotionMetrics_'+m+'/T1_MPR*.txt')
        file = SortFiles(files)[0]
        if d[:2] == 'OF':
            metrics = np.loadtxt(file, unpack=True, usecols=0)
        if d[:2] == 'ON':
            metrics = np.loadtxt(file, unpack=True, usecols=3)

        for name in names:
            s = np.where(subj==name)[0]
            RMS_ref = metrics_ref[s]
            RMS = metrics[s]

            lines.append('Input '+name[:-1]+'.long.X_BASE_'+name[:-1]+' Main '+str(RMS_ref[0]))
            if 'STILL' in d:
                lines.append('Input X_'+d+name[:-1]+'.long.X_BASE_'+name[:-1]+' Main '+str(RMS[0]))
            else:
                lines.append('Input X_'+d+'RR_'+name[:-1]+'.long.X_BASE_'+name[:-1]+' Main '+str(RMS[0]))
        print(d + ' done')

        save = np.array(lines)
        head = 'GroupDescriptorFile 1 \nMeasurementName thickness \nClass Main'
        head += ' \nVariables RMS'
        np.savetxt(outDir+'MOCO_'+d+'long.fsgd', save, header=head, comments='', fmt='%s')



    '''Assemble the data: '''
    # 10 to 15 minutes
    descr = ['MOCO_ON_STILL_long', 'MOCO_ON_NOD_long',  'MOCO_OFF_NOD_long', 'MOCO_ON_SHAKE_long', 'MOCO_OFF_SHAKE_long']
    for d in descr:
        # lh
        # resample to fsaverage and put all in one file:
        subprocess.run('mris_preproc --fsgd '+outDir+d+'.fsgd --target fsaverage --hemi lh --meas thickness --out '+outDir+'lh.'+d+'.thickness.00.mgh', shell=True)
        # smooth each subject's resampled data by 10mm FWHM
        subprocess.run('mri_surf2surf --hemi lh --s fsaverage --sval '+outDir+'lh.'+d+'.thickness.00.mgh --fwhm 10 --cortex --tval '+outDir+'lh.'+d+'.thickness.10.mgh', shell=True)

        # rh
        # resample to fsaverage and put all in one file:
        subprocess.run('mris_preproc --fsgd '+outDir+d+'.fsgd --target fsaverage --hemi rh --meas thickness --out '+outDir+'rh.'+d+'.thickness.00.mgh', shell=True)
        # smooth each subject's resampled data by 10mm FWHM
        subprocess.run('mri_surf2surf --hemi rh --s fsaverage --sval '+outDir+'rh.'+d+'.thickness.00.mgh --fwhm 10 --cortex --tval '+outDir+'rh.'+d+'.thickness.10.mgh', shell=True)



    '''GLM Analysis: '''
    for d in descr:
        print(d)
        # fit glm lh
        subprocess.run('mri_glmfit --y '+outDir+'lh.'+d+'.thickness.10.mgh --fsgd '+outDir+d+'.fsgd --C '+outDir+'lh-Avg-thickness-RMS-Cor.mtx --surf fsaverage lh --cortex --glmdir '+outDir+'lh.'+d+'.glmdir', shell=True)

        # fit glm rh
        subprocess.run('mri_glmfit --y '+outDir+'rh.'+d+'.thickness.10.mgh --fsgd '+outDir+d+'.fsgd --C '+outDir+'rh-Avg-thickness-RMS-Cor.mtx --surf fsaverage rh --cortex --glmdir '+outDir+'rh.'+d+'.glmdir', shell=True)

        # fdr correction (lh and rh together):
        subprocess.run('mri_fdr --i '+outDir+'lh.'+d+'.glmdir/lh-Avg-thickness-RMS-Cor/sig.mgh nomask '+outDir+'lh.'+d+'.glmdir/lh-Avg-thickness-RMS-Cor/sig_fdr.mgh --i '+outDir+'rh.'+d+'.glmdir/rh-Avg-thickness-RMS-Cor/sig.mgh nomask '+outDir+'rh.'+d+'.glmdir/rh-Avg-thickness-RMS-Cor/sig_fdr.mgh --fdr 0.05', shell=True)


        # view the results (significance maps and beta - change)- opens new freeview window for each!!!
        folder = outDir+'lh.'+d+'.glmdir/lh-Avg-thickness-RMS-Cor/'
        subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'sig_fdr.mgh:overlay_threshold=1.301,6 -viewport 3d', shell=True)
        folder = outDir+'lh.'+d+'.glmdir/'
        subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'beta.mgh:overlay_threshold=0.04,0.2 -viewport 3d', shell=True)


        folder = outDir+'rh.'+d+'.glmdir/rh-Avg-thickness-RMS-Cor/'
        subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/rh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'sig_fdr.mgh:overlay_threshold=1.301,6 -viewport 3d', shell=True)
        folder = outDir+'rh.'+d+'.glmdir/'
        subprocess.run('freeview -f $SUBJECTS_DIR/fsaverage/surf/rh.inflated:annot=aparc.annot:annot_outline=1:overlay='+folder+'beta.mgh:overlay_threshold=0.04,0.2 -viewport 3d', shell=True)





