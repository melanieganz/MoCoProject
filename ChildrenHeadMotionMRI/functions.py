'''
Utitlity funcitons for analsing children's motion patterns with the Jupyter 
noteboks 'Matrix_Extraction.ipynb' and 'Motion_Patterns.ipynb'
'''

import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import scipy.io
import glob
from transforms3d.affines import decompose
from transforms3d.euler import mat2euler
import pydicom
import datetime
import csv


def read_nii(path):
    img = nib.load(path)
    return img

# compute centroids based on coordinates of voxels making up the segment of the region
def compute_centroid(indeces):
    cx = np.average(indeces[0])
    cy = np.average(indeces[1])
    cz = np.average(indeces[2])
    return int(cx), int(cy), int(cz)


# path - to the regions binarization
# matrix - vox2ras from freesurfer
# xOffset - used to move centroids left-right
def getCentriods(path,matrix,xOffset=0):
    result = np.zeros((256,256,256))
    result2 = np.zeros((256,256,256))
    print(path)

    minDist = 1000
    centroids = []
    regions = []
    # for region classification counter
    xx = 2
    img1 = 0
    for root, dirs, files in os.walk(path):
        print(files)
        # interate through the 16 binarized files
        for file in files:
            if "region" in file:
                minDist = 1000
                path = root+"/"+file
                regions.append(file)
                img1 = read_nii(path)
                img = img1.get_fdata()
                indeces = np.where(img == 1)
                x,y,z = compute_centroid(indeces)
                xx+=1
                # check if centroid is inside region
                if img[x][y][z]!=1:
                    for a in range(len(indeces[0])):
                        # compute the ssd between obtained centroid and all the region voxels
                        ssd = np.sqrt(((np.array([indeces[0][a],indeces[1][a],indeces[2][a]])-np.array([x,y,z]))**2).sum())
                        if ssd < minDist:
                            minDist = ssd
                            # update centroid with the closest voxel so far 
                            cx,cy,cz = indeces[0][a],indeces[1][a],indeces[2][a]
                    print(x,y,z,cx,cy,cz)        
                else:
                    cx,cy,cz = x,y,z
                point = (np.array([cx,cy,cz,1]))
                # apply vox2ras
                newp = matrix@point 
                centroids.append([newp[0]+xOffset,newp[1],newp[2]])
                if xx<11:
                    for q in range(cx-5,cx+5):
                        for w in range(cy-5,cy+5):
                            for e in range(cz-5,cz+5):
                                result[q][w][e]=1   
                else:
                    for q in range(cx-5,cx+5):
                        for w in range(cy-5,cy+5):
                            for e in range(cz-5,cz+5): 
                                result2[q][w][e]=1   
                    
    new_img = nib.Nifti1Image(result, img1.affine, img1.header)
    new_img.to_filename("subcortical.nii")
    new_img1 = nib.Nifti1Image(result2, img1.affine, img1.header)
    new_img1.to_filename("cortical.nii")
    return centroids,regions,result


def ssd(A,B):
    return np.sqrt(((np.array(A)-np.array(B))**2).sum())


def computeHistograms2(filenames,centroids):
    histograms = []
    for centroid in centroids:
        histogram = []
        for file in filenames:
            mat = scipy.io.loadmat(file)
            transf = mat["matrix"]
            ec = np.append(centroid,1)
            newCentroid = np.dot(transf, ec)
            ssd1 = ssd(ec,newCentroid)
            histogram.append(ssd1)
        histograms.append(histogram)    
    return histograms
    

def shortRegName(regionsMel):
        l = regionsMel.split("(")[-1].split(")")[0].split(" ")
        while "" in l: l.remove("")
        y = l[-1].replace('ctx-','')
        return y 

        
def MakeSubPlot(data1,data2,data3,label1,label2,label3,title,ymax,j,k,nr,cort=0):
    data = []
    for i in range(j,k):
        data.append(data1[:,i])
        data.append(data2[:,i])
        data.append(data3[:,i])
    dataa = np.array(data)
    vox2ras =[[-1.00000,0.00000,0.00000,133.17697],[0.00000,0.00000,1.00000,-105.64328],[0.00000,-1.00000, 0.00000,88.49352],
              [0.00000,0.00000,0.00000,1.00000]]
    centroidsMel, regionsMel,result = getCentriods("/pc_disk1/moco/StudentProjects/MSc/Hannah/FrontRadio_Eichhorn_2021/analysis/test data",vox2ras,0)
    plt.subplot(4,1,nr)
    plt.boxplot(data)
    for i in range((k-j)*3):
        if i<3:
            if i%3==0:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:blue',label = label1)
            elif i%3==1:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:orange',label= label2)
            else:    
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:green',label= label3)

        else:
            if i%3==0:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:blue')
            elif i%3==1:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:orange')
            else:    
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:green')
    x = []

    for i in range(j*3,k*3):
        if i%3==0:
            x.append(shortRegName(regionsMel[i//3]))
        else:
            x.append("")
    plt.xticks(np.arange((k-j)*3)+2,x)  #+1
    if nr == 4:
        plt.legend(loc="upper center", ncol = 3, bbox_to_anchor=(0.5, -0.3))
    plt.title(title,fontsize=15)
    plt.ylim([0,ymax])
    plt.xlabel("Regions",fontsize=15)
    plt.ylabel("Displacement (mm)",fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='minor', labelsize=15)


def makePlotCh(data1,data2,label1,label2,title,ymax,cort=0, lab_mf=False, leg=False):
    data = []
    for i in range(len(data1[0])):
        data.append(data1[:,i])
        data.append(data2[:,i])
    dataa = np.array(data)
    vox2ras =[[-1.00000,0.00000,0.00000,133.17697],[0.00000,0.00000,1.00000,-105.64328],[0.00000,-1.00000, 0.00000,88.49352],
              [0.00000,0.00000,0.00000,1.00000]]
    centroidsMel, regionsMel,result = getCentriods("/pc_disk1/moco/StudentProjects/MSc/Hannah/FrontRadio_Eichhorn_2021/analysis/test data",vox2ras,0)

    plt.boxplot(dataa)
    for i in range(len(data1[0])*2):
        if i<2:
            if i%2==0:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:blue',label = label1)
            else:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:orange',label= label2)
        else:
            if i%2==0:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:blue')
            else:
                plt.scatter(np.full(len(dataa[i]), i+1),dataa[i],color='tab:orange')
    x = []
    if cort:
        for i in range(16):
            if i%2==0:
                x.append(shortRegName(regionsMel[i//2]))
            else:
                x.append("")
    else:
         for i in range(16,32):
            if i%2==0:
                x.append(shortRegName(regionsMel[i//2]))
            else:
                x.append("")
    plt.xticks(np.arange(len(data1[0])*2)+1,x)
    if leg==True:
        plt.legend(bbox_to_anchor=(1.02,1.3), loc='upper left', borderaxespad=0.)
    plt.title(title,fontsize=15)
    plt.ylim([0,ymax])
    plt.xlabel("Regions",fontsize=15)
    plt.ylabel("Displacement (mm)",fontsize=15)
    if lab_mf==True:
        plt.ylabel("Motion-free time (%)",fontsize=15) 
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.tick_params(axis='both', which='minor', labelsize=9)
    
    
    
def search_string_in_file(file_name, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
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
    # Return list of tuples containing line numbers and lines where string is found
    return list_of_results


def PClCentroid(sub):
    track_dir = '/mnt/mocodata1/BUF/DataTCL/'
    track_file = glob.glob(track_dir+sub+'/*MOT.tsm')[0]
    
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
    coord =  np.array([float(x),float(y),float(z),1])
    
    # move into RAS space:
    cal_file = glob.glob(track_dir+sub+'/*ALN.tsa')[0]
    find = search_string_in_file(cal_file, 'A_tcs2dcs')
    if len(find)>0:
        find = find[0]
        skip = find[0]
    else:
        skip=0 
    cal_ = np.loadtxt(cal_file, skiprows=skip)
    cal = np.dot(np.array([[-1,0,0,0], [0,1,0,0], [0,0,-1,0], [0,0,0,1]]), cal_, )
    
    coord_RAS = np.matmul(cal, coord)
    
    return coord_RAS[0:3]


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


def DecomposeMat(matrices):
    '''
    Decomposes all elements of matrices into translation and rotation for one subject.

    Parameters
    ----------
    matrices : np.array, size (N, 4, 4)
        Matrices for all time points in RAS.

    Returns
    -------
    mean_transl, mean_rot
        mean values for translation and rotation on each axis.
    std_transl, std_rot
        standard deviations for translation and rotation on each axis

    '''
    
    transl, rot = np.zeros((len(matrices), 3)), np.zeros((len(matrices), 3))
    for i in range(0, len(matrices)):
        transl[i], rot[i], temp = ParametersFromTransf(matrices[i])
        
    mean_transl = np.mean(abs(transl), axis=0)
    std_transl = np.std(abs(transl), axis=0)
    mean_rot = np.mean(abs(rot), axis=0)
    std_rot = np.std(abs(rot), axis=0)
    
    return mean_transl, mean_rot, std_transl, std_rot


def Decompose_All(subjs, path):
    
    mean_tr, mean_rot, std_tr, std_rot = [], [], [], []    
    for sub in subjs:
        Scan_Times = np.loadtxt('/pc_disk1/moco/StudentProjects/MSc/Hannah/FrontRadio_Eichhorn_2021/analysis/Scan_Times.txt', dtype=str, delimiter=', ', skiprows=1)
        mat_files = StartEndFrameNr([sub], Scan_Times, path)
        #mat_files = os.listdir(path+'matrices '+sub+'/')
        matrices = np.zeros((len(mat_files), 4, 4))
        cnt = 0
        for f in mat_files:
            mat = scipy.io.loadmat(f)
            matrices[cnt] = mat["matrix"]
            cnt+=1
        
        mt, mr, st, sr = DecomposeMat(matrices)
        mean_tr.append(mt)
        mean_rot.append(mr)
        std_tr.append(st)
        std_rot.append(sr)
          
    return np.array(mean_tr), np.array(mean_rot), np.array(std_tr), np.array(std_rot)


def DecomposeMat_pm(matrices):
    '''
    Decomposes all elements of matrices into translation and rotation for one subject.

    Parameters
    ----------
    matrices : np.array, size (N, 4, 4)
        Matrices for all time points in RAS.

    Returns
    -------
    mean_transl, mean_rot
        mean values for translation and rotation on each axis.
    std_transl, std_rot
        standard deviations for translation and rotation on each axis

    '''
    
    transl, rot = np.zeros((len(matrices), 3)), np.zeros((len(matrices), 3))
    for i in range(0, len(matrices)):
        transl[i], rot[i], temp = ParametersFromTransf(matrices[i])
        
    mean_transl = np.mean(transl, axis=0)
    std_transl = np.std(transl, axis=0)
    mean_rot = np.mean(rot, axis=0)
    std_rot = np.std(rot, axis=0)
    
    return mean_transl, mean_rot, std_transl, std_rot


def Decompose_All_pm(subjs, path):
    
    mean_tr, mean_rot, std_tr, std_rot = [], [], [], []    
    for sub in subjs:
        Scan_Times = np.loadtxt('/pc_disk1/moco/StudentProjects/MSc/Hannah/FrontRadio_Eichhorn_2021/analysis/Scan_Times.txt', dtype=str, delimiter=', ', skiprows=1)
        mat_files = StartEndFrameNr([sub], Scan_Times, path)
        #mat_files = os.listdir(path+'matrices '+sub+'/')
        matrices = np.zeros((len(mat_files), 4, 4))
        cnt = 0
        for f in mat_files:
            mat = scipy.io.loadmat(f)
            matrices[cnt] = mat["matrix"]
            cnt+=1
        
        mt, mr, st, sr = DecomposeMat_pm(matrices)
        mean_tr.append(mt)
        mean_rot.append(mr)
        std_tr.append(st)
        std_rot.append(sr)
          
    return np.array(mean_tr), np.array(mean_rot), np.array(std_tr), np.array(std_rot)


def readtcs2dcs(path):
    A_tcs2dcs = []
    with open(path ,'r') as mat:
        lines = mat.read().splitlines()[-5:]
        for line in lines:
            l = []
            elem = line.split("  ")
            for e in elem:
                fff = e.strip()
                if fff!="":
                    l.append(float(fff))
            if len(l)!= 0:        
                A_tcs2dcs.append(l)     
    return np.array(A_tcs2dcs)


# given the tim file, the scan start and the time difference between local and remote time, returns the frame number 
# corresponding to the scan start by adding the time difference th local time to obtain the remote time (the scan time is in remote time)
def findFrameNumber(path, time, timeToAdd):
    splitIndex = 0
    with open(path,'r') as f:
        lines = f.read().splitlines()
        for line in lines:
             if "   Camera Time Stamp" in line:
                splitIndex = lines.index(line)
                #print(splitIndex)
                break
        for line in lines[splitIndex+1:]:         
            lineList = line.split("         ")
            systemTime = lineList[-1].strip()
            time_obj = datetime.datetime.strptime(systemTime, '%H:%M:%S.%f')
            remoteTime = (time_obj+ timeToAdd).strftime('%H:%M:%S.%f')
            rt = remoteTime[0:2]+remoteTime[3:5]+remoteTime[6:8] 
            if(time.split(".")[0]==rt):
                #print(time,systemTime)
                return lineList[1].strip().split(" ")[0]
 
  

# finds time difference between local and remote time, adds it to the local time and returns the frame number corresponding to the scan start
def findFrameNumberInLog(path, time):
    splitIndex = 0
    timeDiff = ""
    with open(path,'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if "Time diff" in line:
                l = line.split(": ")
                timeDiff = l[1].strip()
                timeDiffList = timeDiff.split(":")
                a = timeDiffList[-1].split(".")
                timeToAdd = datetime.timedelta(minutes=int(timeDiffList[1]), seconds=int(a[0]), milliseconds=int(a[1]))
                #print(timeToAdd)
                break
        for line in lines:
            if "Time Stamp" in line:
                splitIndex = lines.index(line)
                #print(splitIndex)
                break
        for line in lines[splitIndex+1:]:    
            lineList = line.split("       ")
            systemTime = lineList[-1].strip()
            time_obj = datetime.datetime.strptime(systemTime, '%H:%M:%S.%f')
            remoteTime = (time_obj+ timeToAdd).strftime('%H:%M:%S.%f')
            rt = remoteTime[0:2]+remoteTime[3:5]+remoteTime[6:8] 
            if(time.split(".")[0]==rt):
                print(time,systemTime)
                #print(lineList)
                test = lineList[1].strip().split(" ")[0]
                if test == "":
                    test = lineList[2].strip().split(" ")[0]
                return test
            
            
#finds the time difference between local and remote time in the LOG file 
#(for patients that have a TIM file, this information is at the bottom of the LOG file)

def getTimeDifference(path):
    with open(path,'r') as f:
        lines = f.read().splitlines()
        line = lines[-2]
        timeDiff = line.split("   ") 
        while "" in timeDiff: timeDiff.remove("")
        #print(timeDiff)
        test = timeDiff[2].replace(" ", "")
        if test == 'TimeDiff':
            line = lines[-1]
            timeDiff = line.split("   ") 
            while "" in timeDiff: timeDiff.remove("")
            #print(timeDiff)
            test =  timeDiff[2].replace(" ", "")
        return test
    
    
def getAcqTimeFromDICOM(path):
    ds = pydicom.filereader.dcmread(path)
    return ds.AcquisitionTime


# extract all file paths from the tracking folder
def getFileNames(folder):
    poa = ""
    aln = ""
    tim = ""
    log = ""
    for root, dirs, files in os.walk(folder):
        for file in files:
            if "ALN" in file:
                aln = root + "/"+file
            if "LOG" in file:
                log = root + "/"+file
            if "POA" in file:
                poa = root + "/"+file
            if "TIM" in file:
                tim = root + "/"+file
    return aln, tim, poa, log                        
  
    
# read a cvs file containing the start time
def getTime(file,subject):
    with open(file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0]==subject:
                return row[-1]
            
def transform(path, A_tcs2ecs,startMPRage,savepath):
    i = 0
    refTransf = np.zeros((4,4)) 
    splitIndex = 0
    with open(path ,'r') as mat:
        lines = mat.read().splitlines() 
        for line in lines:
            #identify line where matrics start via the line above
            if "Frame Number" in line:
                splitIndex = lines.index(line)
                print(splitIndex)
                break
        for line in lines[splitIndex+1:-1]:
            t = []
            elem = line.split("    ")
            while "" in elem: elem.remove("")
            if elem[0].strip()==startMPRage:
                for e in elem[1:]:
                    fff = e.strip()
                    if fff!="":
                        t.append(float(fff))
                if(int(t[-1])!=1)  :
                    del t[-1]
                refTransf = np.array(t).reshape((4,4))
                
                #scipy.io.savemat(savepath+"/"+str(len(str(startMPRage)))+"_"+str(startMPRage)+" ref"'.mat', {'matrix': A_tcs2ecs @ (np.linalg.inv(refTransf)) @ refTransf @ (np.linalg.inv(A_tcs2ecs))})
                refIndex = lines[splitIndex+1:-1].index(line)
                break
        # transform all individual matrices to RAS        
        for line in lines[splitIndex+1:-1:10]: 
                t = [];
                elem = line.split("    ")
                while "" in elem: elem.remove("")
                for e in elem[1:]:
                    fff = e.strip()
                    if fff!="":
                        t.append(float(fff))           
                tracolineTransfMat =np.zeros((4,4))
                if(int(t[-1])!=1)  :
                    del t[-1]
                tracolineTransfMat = np.array(t).reshape((4,4))
                A_ref2src_ecs_ref = A_tcs2ecs @ (np.linalg.inv(tracolineTransfMat)) @ refTransf @ (np.linalg.inv(A_tcs2ecs))
                scipy.io.savemat(savepath+"/"+str(len(str(elem[0].strip())))+"_"+str(elem[0].strip())+'.mat', {'matrix': A_ref2src_ecs_ref})
                i+=1                  
          
    return refTransf


def GetFrameNr(time, path):
    aln, tim, poa, log  = getFileNames(path)
    
    # get the frame for the corresponding time either from tim or log file
    if tim!="":
        timeDifference = getTimeDifference(log)
        timeDiffList = timeDifference.split(":")
        a = timeDiffList[-1].split(".")
        # convert to python timedate
        timeToAdd = datetime.timedelta(minutes=int(timeDiffList[1]), seconds=int(a[0]), milliseconds=int(a[1]))
        refFrame = findFrameNumber(tim, time, timeToAdd)
    elif log!="":
        refFrame = findFrameNumberInLog(log, time)
    else:
        print('Neither log nor tim file')
        
    return refFrame


def StartEndFrameNr(subjs, Scan_Times, dir_matr):
    for sub in subjs:
        # look up frame numbers corrsponding to start and end times:
        ind = np.where(Scan_Times[:,0]==sub)[0]
        st = Scan_Times[ind]
        start_times, end_times = st[:,2], st[:,3]
        start_frames, end_frames = [], []
        for start, end in zip(start_times, end_times):
            path = '/mnt/mocodata1/BUF/DataTCL/'+sub+'/'
            test_start = GetFrameNr(start, path)
            print(test_start)
            test_end = GetFrameNr(end, path)
            if test_start is None:
                save_start = False
                timeToAdd = datetime.timedelta(seconds=30)
                start_ = datetime.datetime.strptime(start, '%H%M%S.%f')
                start_ = (start_+ timeToAdd).strftime('%H%M%S.%f')
                trial = GetFrameNr(start_, path)
                if trial is not None:
                    save_start = True 
                    test_start = trial
            else:
                save_start = True
            if test_end is None:
                save_end = False
                timeToAdd = datetime.timedelta(seconds=30)
                end_ = datetime.datetime.strptime(end, '%H%M%S.%f')
                end_ = (end_- timeToAdd).strftime('%H%M%S.%f')
                trial = GetFrameNr(end_, path)
                if trial is not None:
                    save_end = True 
                    test_end = trial
            else:
                save_end = True
            if save_start == True and save_end == True:
                start_frames.append(test_start)
                end_frames.append(test_end)
        
        print(start_frames, end_frames)
        
        start_frames = np.array(start_frames, dtype=int)
        end_frames = np.array(end_frames, dtype=int)
        
        Frame_Nrs = np.arange(start_frames[0], end_frames[0]+1, step=1)
        for i in range(1, len(start_frames)):
            Frame_Nrs = np.append(Frame_Nrs, np.arange(start_frames[i], end_frames[i]+1, step=1))
        
        names = os.listdir(dir_matr+'matricesREF '+sub)
        names_in_seq = []
        for n in names:
            frame_nr = int(n[2:-4])
            if frame_nr in Frame_Nrs:
                names_in_seq.append(dir_matr+'matricesREF '+sub+'/'+n)
                  
        return names_in_seq
    
    
def StartEndFrameNr2(sub, Scan_Times):
    # look up frame numbers corrsponding to start and end times:
    ind = np.where(Scan_Times[:,0]==sub)[0]
    st = Scan_Times[ind]
    start_times, end_times = st[:,2], st[:,3]
    start_frames, end_frames = [], []
    for start, end in zip(start_times, end_times):
        path = '/mnt/mocodata1/BUF/DataTCL/'+sub+'/'
        test_start = GetFrameNr(start, path)
        test_end = GetFrameNr(end, path)
        if test_start is None:
            save_start = False
            timeToAdd = datetime.timedelta(seconds=30)
            start_ = datetime.datetime.strptime(start, '%H%M%S.%f')
            start_ = (start_+ timeToAdd).strftime('%H%M%S.%f')
            trial = GetFrameNr(start_, path)
            if trial is not None:
                save_start = True 
                test_start = trial
        else:
            save_start = True
        if test_end is None:
            save_end = False
            timeToAdd = datetime.timedelta(seconds=30)
            end_ = datetime.datetime.strptime(end, '%H%M%S.%f')
            end_ = (end_- timeToAdd).strftime('%H%M%S.%f')
            trial = GetFrameNr(end_, path)
            if trial is not None:
                save_end = True 
                test_end = trial
        else:
            save_end = True
        if save_start == True and save_end == True:
            start_frames.append(test_start)
            end_frames.append(test_end)

    start_frames = np.array(start_frames, dtype=int)
    end_frames = np.array(end_frames, dtype=int)
    
    return start_times, end_times, start_frames, end_frames
    

def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)
