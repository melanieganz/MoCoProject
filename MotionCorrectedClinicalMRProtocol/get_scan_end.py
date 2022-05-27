import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import datetime as dt
import os
import glob
import subprocess
import ismrmrd
import ismrmrd.xsd
from motion_estimates import search_string_in_file, Add_Time


subdir = [] 
for i in range(1,10):
    subdir.append('Subject_0'+str(i)+'/')
for i in range(10,20):
    subdir.append('Subject_'+str(i)+'/')
for i in range(20,23):
    subdir.append('Subject_'+str(i)+'/')
sequs = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_TIRM', 'T2STAR', 'DIFF']


get_scan_start = False
get_scan_end = True
date = '07_30'

h5_folder = '../RawDatah5/'


if get_scan_start:
    cont = True
    for sequ in sequs:
        if os.path.exists('/mnt/mocodata1/Data_Analysis/Motion_Estimates/ScanTimes_'+sequ+'_'+date+'.txt'):
            print('Text file with /data1/hannah/Motion_Estimates/ScanTimes_'+sequ+'_'+date+'.txt already exists! If you still want to continue set cont to True.')
            cont = False
            break
            
        with open('/mnt/mocodata1/Data_Analysis/Motion_Estimates/ScanTimes_'+sequ+'_'+date+'.txt', 'a') as f:
                    f.write('#subject_ID name_of_sequence scan_start\n')
                    
    #cont = True
    if cont:
        for sub in subdir:
            for sequ in sequs:
                # for FLAIR, DIFF and T2STAR not all volunteers available:
                if sequ in ['DIFF', 'T2_FLAIR', 'T2STAR']:
                    tmp = os.listdir('../BIDSdata_defaced/'+sub)
                    tmp_ = ''
                    test = sequ
                    if sequ == 'DIFF':
                        test = 'TRACEW_B0'
                    if test not in tmp_.join(tmp):
                        continue
                
                # generate names for the different motion scans:
                names = ['MOCO_OFF_STILL_*', 'MOCO_ON_STILL_*', 'MOCO_OFF_NOD_*', 
                         'MOCO_ON_NOD_*']
                names = [n+sequ for n in names]
                seq_types = ['STILL_', 'STILL_', 'NOD_', 'NOD_']
                seq_types = [s+sequ for s in seq_types]
                
                if sequ=='T1_MPR':
                    names = ['MOCO_OFF_STILL_*', 'MOCO_ON_STILL_*', 'MOCO_OFF_NOD_*', 
                             'MOCO_ON_NOD_*', 'MOCO_OFF_SHAKE_*', 'MOCO_ON_SHAKE_*']
                    names = [n+sequ for n in names]
                    seq_types = ['STILL_', 'STILL_', 'NOD_', 'NOD_']
                    seq_types = [s+sequ for s in seq_types]
                
                
                for name in names:
                    # search for the raw data file:
                    search = name.lower()
                    parts = search.split('*')
                    
                    if sequ == 'T2_FLAIR':
                        parts = name.split('*')
                        parts[0] = parts[0].lower()
                        parts[1] = parts[1][0:1]+ parts[1][1:].lower()

                    out_file = glob.glob(h5_folder+sub+'/*'+parts[0]+'*'+parts[1]+'*.h5')[0]
                        
                    # find studyTime:
                    dset = ismrmrd.Dataset(out_file[:-3]+'_2.h5', 'dataset', create_if_needed=False)
                    header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
                    time = header.studyInformation.studyTime
                    time = time.strftime('%H%M%S.%f')
                    
                    # save studyTime:
                    with open('/mnt/mocodata1/Data_Analysis/Motion_Estimates/ScanTimes_'+sequ+'_'+date+'.txt', 'a') as f:
                        f.write(sub+' '+name+' '+time+'\n')


# Now extract end of acquisition:
ScanTimes = {'STILL_T1_MPR':np.array([4,40]), 'NOD_T1_MPR':np.array([5,12]), 
             'SHAKE_T1_MPR':np.array([5,12]), 'STILL_T1_TIRM':np.array([3,10]), 
             'NOD_T1_TIRM':np.array([3,51]), 'STILL_T2_TSE':np.array([2,30]), 
             'NOD_T2_TSE':np.array([3,6]),
             'STILL_T2_FLAIR':np.array([4,12]), 'NOD_T2_FLAIR':np.array([4,47]),
             'STILL_T2STAR':np.array([2,25]), 'NOD_T2STAR':np.array([2,25]),
             'STILL_EPI_SWI':np.array([0,52]), 'NOD_EPI_SWI':np.array([0,52]),
             'STILL_DIFF':np.array([0,42]), 'NOD_DIFF':np.array([0,42])}    

if get_scan_end:
    cont = True
    for sequ in sequs:
        if os.path.exists('/mnt/mocodata1/Data_Analysis/Motion_Estimates/ScanEndTimes_'+sequ+'_'+date+'.txt'):
            print('Text file with /data1/hannah/Motion_Estimates/ScanEndTimes_'+sequ+'_'+date+'.txt already exists! If you still want to continue set cont to True.')
            cont = False
            break
            
        with open('/mnt/mocodata1/Data_Analysis/Motion_Estimates/ScanEndTimes_'+sequ+'_'+date+'.txt', 'a') as f:
                    f.write('#subject_ID name_of_sequence scan_end\n')
                    
    #cont = True
    if cont:

        for sub in subdir:
            for sequ in sequs:
                # for FLAIR, DIFF and T2STAR not all volunteers available:
                if sequ in ['DIFF', 'T2_FLAIR'] and sub in ['HC_0'+str(i)+'/' for i in range(1,10)]:
                        continue
                if sequ in ['DIFF', 'T2_FLAIR'] and sub in ['HC_'+str(i)+'/' for i in range(10,13)]:
                    continue
                if sequ == 'T2STAR' and sub in ['HC_20/', 'HC_21/', 'HC_22/']:
                    continue
                
                # generate names for the different motion scans:
                names = ['MOCO_OFF_STILL_*', 'MOCO_ON_STILL_*', 'MOCO_OFF_NOD_*', 
                         'MOCO_ON_NOD_*']
                names = [n+sequ for n in names]
                seq_types = ['STILL_', 'STILL_', 'NOD_', 'NOD_']
                seq_types = [s+sequ for s in seq_types]
                
                if sequ=='T1_MPR':
                    names = ['MOCO_OFF_STILL_*', 'MOCO_ON_STILL_*', 'MOCO_OFF_NOD_*', 
                             'MOCO_ON_NOD_*', 'MOCO_OFF_SHAKE_*', 'MOCO_ON_SHAKE_*']
                    names = [n+sequ for n in names]
                    seq_types = ['STILL_', 'STILL_', 'NOD_', 'NOD_', 'SHAKE_', 'SHAKE_']
                    seq_types = [s+sequ for s in seq_types]
                
                
                for name, seq_type in zip(names, seq_types):  
                    file_scan_start = '/mnt/mocodata1/Data_Analysis/Motion_Estimates/ScanTimes_'+sequ+'_'+date+'.txt'
                    
                    # search for sub and name:
                    find = int(search_string_in_file(file_scan_start, sub)[0][0])
                    with open(file_scan_start, 'r') as read_obj:
                        lines = read_obj.readlines()
                    
                    time = None
                    for l in lines[find-1:]:
                        if name in l and sub in l:
                            time = l.split()[2]
                            break
                    if time == None:
                        print(name[:-1]+' cannot be found for '+sub+' in file: '+out_file)
                        print('Does the name conatin an *?')
                            
                    duration = ScanTimes[seq_type]
                    end_time = Add_Time(time, add_min=int(duration[0]), add_sec=int(duration[1]))
                    
                    # save the time:
                    with open('/mnt/mocodata1/Data_Analysis/Motion_Estimates/ScanEndTimes_'+sequ+'_'+date+'.txt', 'a') as f:
                        f.write(sub+' '+name+' '+end_time+'\n')
        

            
            

