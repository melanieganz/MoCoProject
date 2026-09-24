'''
Extracts scan start times (from the raw ismrmrd headers) and scan end
times (start + a fixed per-sequence duration) for every acquisition.

This is entirely optional: the OpenNeuro release (ds004332) already
ships the extracted scan end times under derivatives/results (see the
top-level README), so this script only needs to be run again if those
need to be regenerated from scratch.

Unlike the rest of the pipeline, this script needs the *raw* ismrmrd
data, which is not part of the OpenNeuro/BIDS release. It is available,
on request, from PublicnEUro:
    PN000009 "Markerless Prospective Motion Correction"
    https://datacatalog.publicneuro.eu/dataset/PN000009%20Markerless%20Prospective%20Motion%20Correction/V1

See the README.md in this folder for the required environment
variables and dependencies, and for an important caveat about the raw
subject-numbering assumption this script makes.
'''
import os
import sys
import glob
import numpy as np
import ismrmrd
import ismrmrd.xsd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from motion_estimates import search_string_in_file, Add_Time  # noqa: E402

root = os.environ.get("MOCO_DATASET_PATH")
raw_root = os.environ.get("PUBLICNEURO_RAW_PATH")

out_dir = os.path.join(root, 'derivatives/results/Motion_Estimates/')
os.makedirs(out_dir, exist_ok=True)


def raw_to_bids_subject(raw_subject_dir):
    '''
    Maps a raw PublicnEUro subject directory name (e.g. 'Subject_01/')
    to the corresponding OpenNeuro/BIDS subject ID (e.g. 'sub-01/').

    ASSUMPTION: raw subject numbering matches the OpenNeuro ds004332
    numbering 1:1 in the same order. This has not been verified against
    an actual PublicnEUro download (access is granted upon request) --
    please confirm this mapping is correct before trusting results from
    this script, and update this function if it isn't.
    '''
    num = ''.join(ch for ch in raw_subject_dir if ch.isdigit())
    return f'sub-{int(num):02d}/'


subdir = []
for i in range(1, 10):
    subdir.append('Subject_0'+str(i)+'/')
for i in range(10, 20):
    subdir.append('Subject_'+str(i)+'/')
for i in range(20, 23):
    subdir.append('Subject_'+str(i)+'/')
sequs = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_TIRM', 'T2STAR', 'DIFF']


get_scan_start = False
get_scan_end = True
date = '07_30'


def has_sequence(bids_sub, sequ):
    '''Whether this subject's BIDS anat/ directory has a sidecar for sequ.'''
    tag = 'TRACEW_B0' if sequ == 'DIFF' else sequ
    tmp = os.listdir(os.path.join(root, bids_sub, 'anat'))
    return tag in ''.join(tmp)


if get_scan_start:
    cont = True
    for sequ in sequs:
        if os.path.exists(out_dir + 'ScanTimes_'+sequ+'_'+date+'.txt'):
            print('Text file with '+out_dir+'ScanTimes_'+sequ+'_'+date+'.txt already exists! If you still want to continue set cont to True.')
            cont = False
            break

        with open(out_dir + 'ScanTimes_'+sequ+'_'+date+'.txt', 'a') as f:
                    f.write('#subject_ID name_of_sequence scan_start\n')

    #cont = True
    if cont:
        for sub in subdir:
            bids_sub = raw_to_bids_subject(sub)
            for sequ in sequs:
                # for FLAIR, DIFF and T2STAR not all volunteers available:
                if sequ in ['DIFF', 'T2_FLAIR', 'T2STAR']:
                    if not has_sequence(bids_sub, sequ):
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

                    out_file = glob.glob(raw_root+sub+'/*'+parts[0]+'*'+parts[1]+'*.h5')[0]

                    # find studyTime:
                    dset = ismrmrd.Dataset(out_file[:-3]+'_2.h5', 'dataset', create_if_needed=False)
                    header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
                    time = header.studyInformation.studyTime
                    time = time.strftime('%H%M%S.%f')

                    # save studyTime:
                    with open(out_dir + 'ScanTimes_'+sequ+'_'+date+'.txt', 'a') as f:
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
        if os.path.exists(out_dir + 'ScanEndTimes_'+sequ+'_'+date+'.txt'):
            print('Text file with '+out_dir+'ScanEndTimes_'+sequ+'_'+date+'.txt already exists! If you still want to continue set cont to True.')
            cont = False
            break

        with open(out_dir + 'ScanEndTimes_'+sequ+'_'+date+'.txt', 'a') as f:
                    f.write('#subject_ID name_of_sequence scan_end\n')

    #cont = True
    if cont:

        for sub in subdir:
            bids_sub = raw_to_bids_subject(sub)
            for sequ in sequs:
                # for FLAIR, DIFF and T2STAR not all volunteers available:
                if sequ in ['DIFF', 'T2_FLAIR', 'T2STAR'] and not has_sequence(bids_sub, sequ):
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
                    file_scan_start = out_dir + 'ScanTimes_'+sequ+'_'+date+'.txt'

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
                        print(name[:-1]+' cannot be found for '+sub+' in file: '+file_scan_start)
                        print('Does the name conatin an *?')
                        continue

                    duration = ScanTimes[seq_type]
                    end_time = Add_Time(time, add_min=int(duration[0]), add_sec=int(duration[1]))

                    # save the time:
                    with open(out_dir + 'ScanEndTimes_'+sequ+'_'+date+'.txt', 'a') as f:
                        f.write(sub+' '+name+' '+end_time+'\n')
