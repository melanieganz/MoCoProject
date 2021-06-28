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



def compute_centroid(indeces):
    '''
    Compute centorids based on coordinates of voxels corresponding to the region

    Parameters
    ----------
    indeces : array
        3xN Array of indices corresonding to the region.

    Returns
    -------
    cx, cy, cz : int
        coordinates of the region's centroid.
        
    '''
    cx = np.average(indeces[0])
    cy = np.average(indeces[1])
    cz = np.average(indeces[2])
    return int(cx), int(cy), int(cz)


def getCentroids(path,matrix,xOffset=0):
    '''
    Loads the regions binarizations and calculates a centroid for each region.

    Parameters
    ----------
    path : str
        Path to the regions' binarization.
    matrix : array
        vox2ras matrix from FreeSurfer.
    xOffset : int, optional
        Offset used to move centroids left - right. The default is 0.

    Returns
    -------
    centroids : list
        List of centroid coordinates.
    regions : list
        List of filenames corresponding to each region.

    '''
    result = np.zeros((256,256,256))
    result2 = np.zeros((256,256,256))

    minDist = 1000
    centroids = []
    regions = []
    # for region classification counter
    xx = 2
    for root, dirs, files in os.walk(path):
        # interate through the 16 binarized files
        for file in files:
            if "region" in file:
                minDist = 1000
                path = root+"/"+file
                regions.append(file)
                img1 = nib.load(path)
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
    return centroids,regions


def ssd(A,B):
    '''
    Calculates displacement between arrays A and B.

    Parameters
    ----------
    A : array
        Array.
    B : array
        Array of same shape as A.

    Returns
    -------
    int
        value of sum of squared differences.

    '''
    return np.sqrt(((np.array(A)-np.array(B))**2).sum())


def computeHistograms2(filenames,centroids):
    '''
    Computes displacement of centroids (using ssd) for each matrix from the 
    files specified under filenames.

    Parameters
    ----------
    filenames : list of strings
        List of filenames containing matrices.
    centroids : list
        List of centroid coordinates.

    Returns
    -------
    histograms : list
        List of lists of displacements.

    '''
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
    '''
    Abbreviate the names of regions.

    Parameters
    ----------
    regionsMel : str
        String with region description.

    Returns
    -------
    y : str
        Abbreviation.

    '''
    l = regionsMel.split("(")[-1].split(")")[0].split(" ")
    while "" in l: l.remove("")
    y = l[-1].replace('ctx-','')
    return y 

        
def MakeSubPlot(data1,data2,data3,label1,label2,label3,title,ymax,j,k,nr,cort=0):
    '''
    Automisation of plotting three different metrics for each region.

    Parameters
    ----------
    data1 : array
        First metric.
    data2 : array
        Second metric.
    data3 : array
        Thrid metric.
    label1 : str
        Label of first metric.
    label2 : str
        Label of second metric.
    label3 : str
        label of third metric.
    title : str
        title of subplot.
    ymax : float
        limit for y axis.
    j : int
        select from which region to start iteration.
    k : int
        select until which region to iterate.
    nr : int
        number of subplot.
    cort : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    None.

    '''
    data = []
    for i in range(j,k):
        data.append(data1[:,i])
        data.append(data2[:,i])
        data.append(data3[:,i])
    dataa = np.array(data)
    vox2ras =[[-1.00000,0.00000,0.00000,133.17697],[0.00000,0.00000,1.00000,-105.64328],[0.00000,-1.00000, 0.00000,88.49352],
              [0.00000,0.00000,0.00000,1.00000]]
    centroidsMel, regionsMel = getCentroids("./test data/mel",vox2ras,0)
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
    plt.xticks(np.arange((k-j)*3)+1,x)
    if nr == 4:
        plt.legend(loc="upper center", ncol = 3, bbox_to_anchor=(0.5, -0.3))
    plt.title(title,fontsize=15)
    plt.ylim([0,ymax])
    plt.xlabel("Regions",fontsize=15)
    plt.ylabel("Displacement (mm)",fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='minor', labelsize=15)


def makePlotCh(data1,data2,label1,label2,title,ymax,cort=0, lab_mf=False, leg=False):
    '''
    Plot two datasets comparatively as boxplots.

    Parameters
    ----------
    data1 : array
        first dataset.
    data2 : array
        second dataset.
    label1 : str
        label for first dataset.
    label2 : str
        label for second dataset.
    title : str
        label for title of theplot.
    ymax : float
        limit for y axis.
    cort : int, optional
        Defines whether cortical or subcortical region. The default is 0.
    lab_mf : bool, optional
        Whether y label is displacement of motion free time. The default is False.
    leg : bool, optional
        If legend should be added. The default is False.

    Returns
    -------
    None.

    '''
    data = []
    for i in range(len(data1[0])):
        data.append(data1[:,i])
        data.append(data2[:,i])
    dataa = np.array(data)
    vox2ras =[[-1.00000,0.00000,0.00000,133.17697],[0.00000,0.00000,1.00000,-105.64328],[0.00000,-1.00000, 0.00000,88.49352],
              [0.00000,0.00000,0.00000,1.00000]]
    centroidsMel, regionsMel = getCentroids("./test data/mel",vox2ras,0)

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
    '''
    Search for the given string in file and return lines containing that string,
    along with line numbers

    Parameters
    ----------
    file_name : str
        name of the file to search.
    string_to_search : str
        string to search in the file.

    Returns
    -------
    list_of_results : list
        list of all ines containing the string.

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
    # Return list of tuples containing line numbers and lines where string is found
    return list_of_results


def PClCentroid(sub):
    '''
    Extract position of Point Cloud Centroid in RAS coordinates from tracking
    files.

    Parameters
    ----------
    sub : str
        ID for the corresponding subject.

    Returns
    -------
    array
        coordinates of point cloud centroid.

    '''
    track_dir = '/mnt/mocodata1/MoCoOther/Andreea-Project/moco/subjects data/'
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
    Decomposes all elements of matrices into translation and rotation for 
    one subject. Averages absolute values.

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
    '''
    Decompose all matrices for all subjects specified in subjs. Averages
    absolute values.

    Parameters
    ----------
    subjs : list of str
        List of all subject IDs for which the matrices should be decomposed.
    path : str
        path to the files containing matrices.

    Returns
    -------
    mean_transl, mean_rot
        list of mean values for translation and rotation on each axis.
    std_transl, std_rot
        list of standard deviations for translation and rotation on each axis

    '''
    
    mean_tr, mean_rot, std_tr, std_rot = [], [], [], []    
    for sub in subjs:
        Scan_Times = np.loadtxt('./new_calc/Scan_Times.txt', dtype=str, delimiter=', ', skiprows=1)
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
    Decomposes all elements of matrices into translation and rotation for 
    one subject. Averages positive and negative values-

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
    '''
    Decompose all matrices for all subjects specified in subjs. Averages
    positive and negative values.

    Parameters
    ----------
    subjs : list of str
        List of all subject IDs for which the matrices should be decomposed.
    path : str
        path to the files containing matrices.

    Returns
    -------
    mean_transl, mean_rot
        list of mean values for translation and rotation on each axis.
    std_transl, std_rot
        list of standard deviations for translation and rotation on each axis

    '''
    
    mean_tr, mean_rot, std_tr, std_rot = [], [], [], []    
    for sub in subjs:
        Scan_Times = np.loadtxt('./new_calc/Scan_Times.txt', dtype=str, delimiter=', ', skiprows=1)
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
    '''
    Loads cross calibration matrix between tracking and scanner coordinate
    system.

    Parameters
    ----------
    path : str
        path to the corresponding file.

    Returns
    -------
    array
        calibration transform.

    '''
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


def findFrameNumber(path, time, timeToAdd):
    '''
    Returns the frame number corresponding to the scan start by adding the 
    time difference to the local time to obtain the remote time

    Parameters
    ----------
    path : str
        path to tim file.
    time : str
        time of scan start.
    timeToAdd : datetime timedelta
        time difference between local and remote time.

    Returns
    -------
    int
        frame number of scan start.

    '''
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
 
  

def findFrameNumberInLog(path, time):
    '''
    Finds the time difference between local and remote time, adds it to the 
    local time and returns the frame number corresponding to the scan start

    Parameters
    ----------
    path : str
        path to log file.
    time : str
        time of scan start.

    Returns
    -------
    int
        frame number of scan start.

    '''
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
            

def getTimeDifference(path):
    '''
    Finds the time difference between local and remote time in the log file,
    for patients with a TIM file, this information is at the bottom of the 
    log file

    Parameters
    ----------
    path : str
        path to log file.

    Returns
    -------
    str
        time difference.

    '''
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
    '''
    Reads the tag corresponding to acquisition time from the DICOM header

    Parameters
    ----------
    path : str
        path to DICOM.

    Returns
    -------
    str
        acquisition time.

    '''
    ds = pydicom.filereader.dcmread(path)
    return ds.AcquisitionTime


def getFileNames(folder):
    '''
    Extracts all relevant file paths (.aln, .tim, .poa, .log) from the folder
    containting tracking data.

    Parameters
    ----------
    folder : str
        path to the folder.

    Returns
    -------
    aln : str
        filename for .aln file.
    tim : str
        filename for .tim file.
    poa : str
        filename for .poa file.
    log : str
        filename for .log file.

    '''
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
  
    
def getTime(file,subject):
    '''
    Reads the custom-made csv file containing the start times and extracts 
    the relevant start time

    Parameters
    ----------
    file : str
        filename of csv file.
    subject : str
        subject ID.

    Returns
    -------
    str
        start time.

    '''
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
    '''
    Extracts the frame number corresponding to the start of scan either from
    .tim or from .log file

    Parameters
    ----------
    time : str
        time of scan start.
    path : str
        path to folder with tracking files.

    Returns
    -------
    refFrame : int
        frame number corresponding to start of scan.

    '''
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
    '''
    Extracts all matrix filenames during MR acquisition (after scan start and 
    before scan end) for all subjects

    Parameters
    ----------
    subjs : list
        list of subject IDs.
    Scan_Times : array
        array of scan times for each subject.
    dir_matr : str
        path to directory with matrices for all subjects.

    Returns
    -------
    names_in_seq : list
        List of all matrix filenames.

    '''
    for sub in subjs:
        # look up frame numbers corrsponding to start and end times:
        ind = np.where(Scan_Times[:,0]==sub)[0]
        st = Scan_Times[ind]
        start_times, end_times = st[:,2], st[:,3]
        start_frames, end_frames = [], []
        for start, end in zip(start_times, end_times):
            path = '/mnt/mocodata1/MoCoOther/Andreea-Project/moco/subjects data/'+sub+'/'
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
    '''
    Extracts start and end times as well as start and end frame numbers for 
    the subjects sub

    Parameters
    ----------
    sub : str
        subject ID.
    Scan_Times : array
        array of scan times for each subject.

    Returns
    -------
    start_times : str
        time of scan start.
    end_times : str
        time of scan end.
    start_frames : int
        frame number of scan start.
    end_frames : int
        frame number of scan end.

    '''
    # look up frame numbers corrsponding to start and end times:
    ind = np.where(Scan_Times[:,0]==sub)[0]
    st = Scan_Times[ind]
    start_times, end_times = st[:,2], st[:,3]
    start_frames, end_frames = [], []
    for start, end in zip(start_times, end_times):
        path = '/mnt/mocodata1/MoCoOther/Andreea-Project/moco/subjects data/'+sub+'/'
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
    '''
    Annotate barplot with p-values

    Parameters
    ----------
    num1 : int
        number of left bar to put bracket over.
    num2 : int
        number of right bar to put bracket over.
    data : str or float
        string to write or number for generating asterixes.
    center : array
        centers of all bars (like plt.bar() input).
    height : array
        heights of all bars (like plt.bar() input).
    yerr : array, optional
        yerrs of all bars (like plt.bar() input). The default is None.
    dh : float, optional
        height offset over bar / bar + yerr in axes coordinates (0 to 1). The 
        default is .05.
    barh : float, optional
        bar height in axes coordinates (0 to 1). The default is .05.
    fs : float, optional
        font size. The default is None.
    maxasterix : int, optional
        maximum number of asterixes to write (for very small p-values). The 
        default is None.

    Returns
    -------
    None.

    '''
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