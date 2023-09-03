import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pydicom
import glob
import json
import datetime as dt
from transforms3d.affines import decompose
from transforms3d.euler import mat2euler


# define scan times for different sequences:
ScanTimes = {'STILL_T1_MPR':np.array([4,40]), 'NOD_T1_MPR':np.array([5,12]),       #Ved ikke hvornår der refereres til ScanTimes
             'SHAKE_T1_MPR':np.array([5,12]), 'STILL_T1_TIRM':np.array([3,10]), 
             'NOD_T1_TIRM':np.array([3,51]), 'STILL_T2_TSE':np.array([2,30]), 
             'NOD_T2_TSE':np.array([3,6]),
             'STILL_T2_FLAIR':np.array([4,12]), 'NOD_T2_FLAIR':np.array([4,47]),
             'STILL_T2STAR':np.array([2,25]), 'NOD_T2STAR':np.array([2,25]),
             'STILL_EPI_SWI':np.array([0,52]), 'NOD_EPI_SWI':np.array([0,52])}
             #'STILL_DIFF':np.array([0,42]), 'NOD_DIFF':np.array([0,42])} #ignore



def ParametersFromTransf(A):
    '''
    Use python module transforms3d to extract translation and rotation 
    parameters from transformation matrix.

    Parameters
    ----------
    A : numpy array (4x4)
        transformation matrix.

    Returns
    -------
    T : numpy array
        translation parameters.
    R : numpy array
        rotation angles in degrees.

    '''
    T, R_, Z_, S_ = decompose(A)
    al, be, ga = mat2euler(R_)
    R = np.array([al*180/np.pi, be*180/np.pi, ga*180/np.pi])
    
    return np.array(T), R, R_


def GetAcquFromJSON(filename):
    '''
    Extract the acquisition time from the JSON document.

    Parameters
    ----------
    filename : str
        path to the JSON file.

    Returns
    -------
    str
        acquisition time.

    '''
    
    with open(filename, 'r') as f:
        data = f.read()
    obj = json.loads(data)
    
    return obj['AcquisitionTime']
    

def Add_Time(start_time, add_min=5, add_sec=0):
    '''
    Parameters
    ----------
    start_time : str
        Start time of acquisition in format %H%M%S.%f.
    add_min : int, optional
        Number of minutes to be added/subtracted to/from start_time. The default is 5.
    add_sec : int, optional
        Number of seconds to be added/subtracted to/from start_time. The default is 0.

    Returns
    -------
    str
        Sting, where add_min and add_sec have been added to start_time.
    '''
    
    time_dt = dt.datetime.strptime(start_time, "%H%M%S.%f")
    end_dt = time_dt+dt.timedelta(minutes=add_min, seconds=add_sec)
    if add_sec < 0 and add_min < 0:
        end_dt = time_dt-dt.timedelta(minutes=add_min, seconds=add_sec)
    
    return dt.datetime.strftime(end_dt, "%H%M%S.%f")


def search_string_in_file(file_name, string_to_search):
    '''
    Search for the given string in file and return lines containing that string,
    along with line numbers

    Parameters
    ----------
    file_name : str
        filename of file to be searched in.
    string_to_search : str
        string that should be searched for.

    Returns
    -------
    list_of_results : list
        list of tuples containing line numbers and lines where string is found.

    '''
    line_number = 0
    list_of_results = []
    
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if string_to_search in line:
                # If yes, then add the line number & line as a tuple in the list
                list_of_results.append((line_number, line.rstrip()))

    return list_of_results


def search_two_strings_in_file(file_name, strings_to_search):
    '''
    Search for two different strings in the filename. Only one of the given 
    strings should apear in the file.

    Parameters
    ----------
    file_name : str
        filename of file to be searched in.
    strings_to_search : list
        list of strings to be searched for.

    Returns
    -------
    tmp : list
        list of tuples containing line numbers and lines where string is found.

    '''
    
    for string in strings_to_search:
        tmp = search_string_in_file(file_name, string)
        if len(tmp)>0:
            return tmp
        

def FindALN(subj, name, bids_dir=None, track_dir=None):
    '''
    Finds the most recent cross calibration for a given sequence.

    Parameters
    ----------
    subj : str
        subject ID.
    name : str
        name of sequence.
    bids_dir : str (optional)
        directory of DICOM files. If None, then it is set to default: 
        '/mnt/mocodata1/Data_Analysis/DICOMS/'+subj
    track_dir : str (optional)
        directory of motion tracking files. If None, then it is set to default: 
        '/mnt/mocodata1/Data_Analysis/TCLData/'+subj

    Returns
    -------
    aln_file : str
        filename of ALN file.

    '''
    if bids_dir is None:
        bids_dir = '../BIDSdata_defaced/'+subj  
    if track_dir is None:
        track_dir = '../TCLData/'+subj
    
    all_aln = glob.glob(track_dir+'*ALN*.tsa')
    file = glob.glob(bids_dir+'*'+name+'*.json')[0] 
    acqu_time = GetAcquFromJSON(file)
    
    diff, aln_file = [], []
    for a in all_aln:
        # find out time difference between remote and local time:
        tmp = a.find('ALN')
        tim = a[0:tmp]+'TIM.tst'
        tmp1, tmp2, local, remote = np.loadtxt(tim, skiprows=11, dtype=str, unpack=True)
        if acqu_time > remote[0] and acqu_time < remote[-1]:
            local_acqu = local[remote>acqu_time][0]
        else:
            local_acqu = acqu_time
        local_acqu = local_acqu[0:2]+local_acqu[3:5]+local_acqu[6:]
        
        find = search_string_in_file(a, 'CrossCalDateAndTime')[0][1]
        time = find[-8:]
        time = time[0:2]+time[3:5]+time[6:]
        if (float(time)-float(local_acqu))<0:
            diff.append(np.abs(float(time)-float(local_acqu)))
            aln_file.append(a)
    ind = np.argmin(diff)

    return aln_file[ind]


def FindPOA_TIM(subj, name, bids_dir=None, track_dir=None):
    '''
    Find the current POA, TIM and MOT files corresponding to the given sequence.

   Parameters
    ----------
    subj : str
        subject ID.
    name : str
        name of sequence.
    bids_dir : str (optional)
        directory of DICOM files. If None, then it is set to default: 
        '../BIDSdata_defaced/'+subj
    track_dir : str (optional)
        directory of motion tracking files. If None, then it is set to default: 
        '../TCLData/'+subj

    Returns
    -------
    poa : str
        POA file (matrices).
    tim : str
        TIM file (time and frame number).
    mot : str
        MOT file (sequence labels and start times).

    '''
    
    if bids_dir is None:
        bids_dir = '../BIDSdata_defaced/'+subj   
    if track_dir is None:
        track_dir = '../TCLData/'+subj
    
    all_tim = glob.glob(track_dir+'*TIM*.tst')
    file = glob.glob(bids_dir+'*'+name+'*.json')[0] 
    acqu_time_j = GetAcquFromJSON(file)
    acqu_time = float(acqu_time_j[0:2]+acqu_time_j[3:5]+acqu_time_j[6:])
    #json file for time XX:XX:0Y.XXX is saved as XX:XX:Y.XXX
    if float(acqu_time_j[6:])<10:
        acqu_time = float(acqu_time_j[0:2]+acqu_time_j[3:5]+'0'+acqu_time_j[6:])
    
    
    for t in all_tim:
        find = int(search_string_in_file(t, 'Point Cloud Number')[0][0])
        with open(t, 'r') as read_obj:
            lines = read_obj.readlines()
        tmp1 = lines[find]
        tmp2 = lines[-2]
        start = float(tmp1[-13:][0:2]+tmp1[-13:][3:5]+tmp1[-13:][6:])
        end = float(tmp2[-13:][0:2]+tmp2[-13:][3:5]+tmp2[-13:][6:])
                       
        if acqu_time > start and acqu_time < end:
            tim = t
            
    poa = tim[0:-7]+'POA.tsp'
    mot = tim[0:-7]+'MOT.tsm'

    return poa, tim, mot

    
def GetTimeFromTCL(mot, name, subj=None, track_dir=None):
    '''
    Searches the scan start in the '_MOT.tsm' file. If the label is available,
    it returns the start stime, otherwise it returns None.

    Parameters
    ----------
    mot : str
        filename of the '_MOT.tsm' file.
    name : str
        name describing the sequence of interest.
    subj : str, optional
        subject ID. Only needs to be provided if track_dir is not provided.
    track_dir : str, optional
        Path to the folder with tracking data. The default is None.

    Returns
    -------
    time : str
        either the time of scan start or None, if the label is not available
        in the '_MOT.tsm' file.

    '''
    
    if track_dir is None:
        track_dir = '../BIDSdata_defaced/'+subj  
    
    find = int(search_string_in_file(mot, 'Label Position')[0][0])
    
    with open(mot, 'r') as read_obj:
        lines = read_obj.readlines()
    time = None
    for l in lines[find:]:
        tmp = [0,0]
        if '*' in name:
            parts = name.split('*')
            if 'T2_FLAIR' in name:
                test = parts[1][0:1]+ parts[1][1:].lower()
                if parts[0].lower() in l and test in l:
                    tmp = l.split()  
                    if tmp[1] == '34':
                        time = tmp[0]
            else:    
                if parts[0].lower() in l and parts[1].lower() in l:
                    tmp = l.split()  
                    if tmp[1] == '34':
                        time = tmp[0]
        else:
            if 'T2_FLAIR' in name:
                n = name.find('T2')
                test = name[:n].lower()+name[n:n+1]+name[n+1:].lower()
                if test in l:
                    tmp = l.split()
                    if tmp[1] == '34':
                        time = tmp[0]         
            else:
                if name.lower() in l:
                    tmp = l.split()
                    if tmp[1] == '34':
                        time = tmp[0]
    
    if time is not None:
        time = time[0:2]+time[3:5]+time[6:]
    
    return time


def GetTimeFromTxt(sub, name):
    '''
    Get the scan end for each volunteer and sequence from a text file
    with times that have been extracted from the raw data tag.

    Parameters
    ----------
    sub : str
        subject ID.
    name : str
        name describing the sequence of interest. Needs to contain a *.

    Returns
    -------
    time : str
        start time for the scan specified with sub and name.

    '''
    
    #find the file with the corresponding scan scart:
    sequ = name.split('*')[1]
    file = '../TCLData/ScanEndTimes_'+sequ+'07_30.txt'
    
    # search for sub and name:
    find = int(search_string_in_file(file, sub)[0][0])
    with open(file, 'r') as read_obj:
        lines = read_obj.readlines()
    
    time = None
    for l in lines[find-1:]:
        if name[:-1] in l and sub in l:
            time = l.split()[2]
            break
    if time == None:
        print(name[:-1]+' cannot be found for '+sub+' in file: '+file)
        print('Does the name contain an *?')

    return time


def FindFrameNr(tim, subj, name, seq_type, mot, bids_dir=None, track_dir=None):
    '''
    Find frame number corresponding to start and end of acquisition of the
    sequence specified in name

    Parameters
    ----------
    tim : str
        TIM.tst filename.
    subj : str
        subject ID.
    name : str
        name describing the sequence of interest.
    seq_type : str
        description for look-up table scan times.
    mot : str
        filename of the '_MOT.tsm' file.
    bids_dir : str, optional
        directory of dicom files. The default is None.
    track_dir : str, optional
        Path to the folder with tracking data. The default is None.

    Returns
    -------
    frame_start : str
        frame number corresponding to start of acquisition.
    frame_end : str
        frame number corresopnding to end of acquisition.

    '''
    
    if bids_dir is None:
        bids_dir = '../BIDSdata_defaced/'+subj  
    
    file = glob.glob(bids_dir+'*'+name+'*.json')[0] 
    acqu_time_j = GetAcquFromJSON(file)
    #json file for time XX:XX:0Y.XXX is saved as XX:XX:Y.XXX
    if float(acqu_time_j[6:])<10:
        acqu_time = acqu_time_j[0:6]+'0'+acqu_time_j[6:]
    else:
        acqu_time = acqu_time_j
    
    end_time = GetTimeFromTxt(subj, name)
    
    # add time to start of acquisition:
    start_time = acqu_time
    #start_time = start_time[0:2]+':'+start_time[2:4]+':'+start_time[4:]
    end_time = end_time[0:2]+':'+end_time[2:4]+':'+end_time[4:]
    # print(start_time, end_time)
    
    frames, tmp1, tmp2, remote = np.loadtxt(tim, skiprows=11, dtype=str, unpack=True)
    frame_start = frames[remote>start_time][0]
    frame_end = frames[remote < end_time][-1]
    
    return frame_start, frame_end


def GetMatrices(tim, poa, mot, subj, name, seq_type, bids_dir=None):
    '''
    Import matrices from poa file (only those during acquisition of sequence 
    specified in name)

    Parameters
    ----------
    tim : str
        TIM.tst filename.
    poa : str
        POA.tsp filename.
    mot : str
        MOT.tsm filename.
    subj : str
        subject ID.
    name : str
        name describing the sequence of interest.
    seq_type : str
        description for look-up table scan times.
    bids_dir : str, optional
        directory of dicom files. The default is None.

    Returns
    -------
    frame_nr : array
        frame numbers for sequence specified in name.
    mat : array
        transformation matrices for sequence specified in name.

    '''
    
    frame_start, frame_end = FindFrameNr(tim, subj, name, seq_type, mot, bids_dir)
    frame_start = int(frame_start)
    frame_end = int(frame_end)
    
    find = search_string_in_file(poa, 'Frame Number')[0]
    skip_ = find[0]
    tmp = np.loadtxt(poa, skiprows=skip_, unpack=True)
    frame_nr = tmp[0].astype(int)
    mat_ = tmp[1:17].T
    mat = mat_.reshape(np.shape(mat_)[0], 4, 4)
    
    # only frame numbers between frame_start and frame_end:
    mat = mat[frame_nr>frame_start]
    frame_nr = frame_nr[frame_nr>frame_start]
    mat = mat[frame_nr<frame_end]
    frame_nr = frame_nr[frame_nr<frame_end]
    
    frames, tmp1, tmp2, remote = np.loadtxt(tim, skiprows=11, dtype=str, unpack=True)
    times = []
    frames = frames.astype(int)
    for f in frame_nr:
        ind = np.where(frames==f)[0][0]
        times.append(remote[ind])
    
    return frame_nr, mat, times




def ExtractMotionParForScan(subj, name, seq_type):
    '''
    Extract Motion parameters for a specific scan

    Parameters
    ----------
    subj : string
        name of subject / subject ID.
    name : string
        scan for which to extract the motion data.
    seq_type : string
        type of sequence (e.g. 'STILL_T1_MPR') for determining scan time

    Returns
    -------
    times
        time points of motion data.
    transl
        translation parameters for each time point specified in times
    rot
        rotation parameters for each time point specified in times

    '''
    
    bids_dir = '../BIDSdata_defaced/'+subj   
    track_dir = '../TCLData/'+subj
    
    name_ = name
    if 'DIFF' in name:
        name_ = name.replace('DIFF', 'TRACEW_B0')
    mat_file, time_file, mot_file = FindPOA_TIM(subj, name_, bids_dir, track_dir)
    
    frame_nr, mat, times = GetMatrices(time_file, mat_file, mot_file, subj, 
                                       name, seq_type, bids_dir)
    
    
    
    cal_file = FindALN(subj, name)
    
    cal_ = np.loadtxt(cal_file, skiprows=16)
    cal = np.dot(np.array([[-1,0,0,0], [0,1,0,0], [0,0,-1,0], [0,0,0,1]]), cal_, )
    
    mat_RAS = np.matmul(cal, np.matmul(np.linalg.inv(mat), 
                                       np.matmul(mat[0], np.linalg.inv(cal))))
    
    transl, rot = np.zeros((len(mat_RAS), 3)), np.zeros((len(mat_RAS), 3))
    for i in range(0, len(mat_RAS)):
        transl[i], rot[i], temp = ParametersFromTransf(mat_RAS[i])
    
    return times, transl, rot


def ExtractMotionMatForScan(subj, name, seq_type, track_dir=None, bids_dir=None):
    '''
    Extract Motion matrices and centroid coordinates (rist time point) for a 
    specific scan

    Parameters
    ----------
    subj : string
        name of subject / subject ID.
    name : string
        scan for which to extract the motion data, use 'None' to extract the motion data for the whole session.
    seq_type : string
        type of sequence (e.g. 'STILL_T1_MPR') for determining scan time
    bids_dir : str (optional)
        directory of DICOM files. If None, then it is set to default: 
        '/mnt/mocodata1/Data_Analysis/DICOMS/'+subj
    track_dir : str (optional)
        directory of motion tracking files. If None, then it is set to default: 
        '/mnt/mocodata1/Data_Analysis/TCLData/'+subj

    Returns
    -------
    coord
        Point cloud centroid for first time point.
    times
        time points of motion data
    mat
        transformation matrix for each time point specified in times, not in 
        RAS system!

    '''
    if bids_dir is None:
        bids_dir = '../BIDSdata_defaced/'+subj   
    if track_dir is None:
        track_dir = '../TCLData/'+subj
    
    name_ = name
    if 'DIFF' in name:
        name_ = name.replace('DIFF', 'TRACEW_B0')
    mat_file, time_file, mot_file = FindPOA_TIM(subj, name_, bids_dir, track_dir)
    
    frame_nr, mat, times = GetMatrices(time_file, mat_file, mot_file, subj, 
                                       name, seq_type, bids_dir)
    
    track_file = glob.glob(track_dir+'*MOT.tsm')[0]
    
    # centroid position:
    find = search_string_in_file(track_file, 'Centroid coordinates')[0]
    skip = find[0]-1
    tmp = np.loadtxt(track_file, skiprows=skip, max_rows=1, dtype=str)
    tmp = tmp[2][1:-1]
    ind = np.char.find(tmp, ',')
    x = tmp[0:ind]
    tmp = tmp[ind+1:]
    ind = np.char.find(tmp, ',')
    y = tmp[0:ind]
    z = tmp[ind+1:]
    coord = np.array([float(x),float(y),float(z)])
    
    return coord, times, mat


def CalcMotionMetricsforScan(subj, name, seq_type, track_dir=None, bids_dir=None, center=False, average=False):
    '''
    Loads the motion data (matrices) and calculates the motion metrics for 
    a sequence specified in name and seq_type.

    Parameters
    ----------
    subj : string
        name of subject / subject ID.
    name : string
        scan for which to extract the motion data, use 'None' to extract the motion data for the whole session.
    seq_type : string
        type of sequence (e.g. 'STILL_T1_MPR') for determining scan time
    bids_dir : str (optional)
        directory of DICOM files. If None, then it is set to default: 
        '/mnt/mocodata1/Data_Analysis/DICOMS/'+subj
    track_dir : str (optional)
        directory of motion tracking files. If None, then it is set to default: 
        '/mnt/mocodata1/Data_Analysis/TCLData/'+subj
    center : bool (optional)
        determines whether the displacements are calculated with respect to 
        the start of the sequence (default) or if True, with respect to the 
        center of acquisition

    Returns
    -------
    magn : array
        displacements for each time point.
    RMS : float
        RMS displacement.
    med_disp : float
        median displacement.
    max_disp : float
        maximum displacement.

    '''
    
    coord_, times, mat = ExtractMotionMatForScan(subj=subj, name=name, 
                                                seq_type=seq_type)
    
    # apply the matrices to the centroid coord:
    coord = np.ones((4))
    coord[0:3] = coord_
    tr_coord = np.matmul(mat, coord)
    
    
    #calculate displacement and its magnitude:
    displ = tr_coord - tr_coord[0]   
    # with respect to center:
    if center == True:
        n = int(len(tr_coord)/2)
        displ = tr_coord - tr_coord[n]
    if average == True:
        av = np.mean(tr_coord, axis=0)
        displ = tr_coord - av
    displ = displ[:,:3]
    magn = np.sqrt(displ[:,0]**2+displ[:,1]**2+displ[:,2]**2)
    
    # calculate RMS:
    n = len(magn)
    RMS = np.sqrt(np.sum(magn**2/n))
    
    # calculate median and max for each scan:
    med_disp = np.median(magn)
    max_disp = np.max(magn)
    
    return magn, RMS, med_disp, max_disp


