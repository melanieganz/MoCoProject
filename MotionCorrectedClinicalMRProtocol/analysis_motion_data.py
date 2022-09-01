'''Warning - this script loads a bunch of motion data points and therefore takes a while to run e.g. over night'''
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import glob
import os
from statsmodels.stats.multitest import multipletests
from motion_estimates import CalcMotionMetricsforScan
from utils import Show_Stars, SortFiles, MakeBoxplot
from statistical_tests import PerformWilcoxonMotion



# Define directories, sequences and volunteers to analyse:
out_dir = '../Results/Motion_Estimates/'

subdir = [] 
for i in range(1,10):
    subdir.append('Subject_0'+str(i)+'/')
for i in range(10,20):
    subdir.append('Subject_'+str(i)+'/')
for i in range(20,23):
    subdir.append('Subject_'+str(i)+'/')
sequs = ['T1_MPR', 'T2_TSE', 'T1_TIRM', 'T2_FLAIR', 'T2STAR', 'DIFF']


save = '_2022_06_02'   
new_calc = True
plot = False


''' (1) Calculate the metrics for each sequence and volunteer: '''

if new_calc:
    for sequ in sequs:
        motion = ['STILL', 'NOD']
        if sequ == 'T1_MPR':
            motion = ['STILL', 'NOD', 'SHAKE']
        
        for mot in motion:
            save_list = []
            print(sequ, mot)
            for subj in subdir:
                ##DEBUG
                print('DEBUG:')
                print(subj)
                
                # skip all volunteers where the following sequences were not acquired
                if sequ in ['DIFF', 'T2_FLAIR', 'T2STAR']:
                    #print('inside first if')
                    tmp = os.listdir('../BIDSdata_defaced/'+subj+'/')
                    tmp_ = ''
                    test = sequ
                    if sequ == 'DIFF':
                        test = 'TRACEW_B0'
                    if test not in tmp_.join(tmp):
                        continue
                
                # MoCo On:
                name = 'MOCO_ON_'+mot+'_*'+sequ+'_'   
                seq_type = mot+'_'+sequ    
                if sequ == 'TRACEW_B0':
                    seq_type = mot+'_DIFF'
                            
                magn_On, RMS_On, med_disp_On, max_disp_On = CalcMotionMetricsforScan(subj, name, 
                                                                                     seq_type)
                
                # MoCo Off:
                name = 'MOCO_OFF_'+mot+'_*'+sequ+'_'   
                seq_type = mot+'_'+sequ    
                if sequ == 'TRACEW_B0':
                    seq_type = mot+'_DIFF'
                
                magn_Off, RMS_Off, med_disp_Off, max_disp_Off = CalcMotionMetricsforScan(subj, 
                                                                                         name, 
                                                                                         seq_type)
                
                save_list.append([RMS_Off, med_disp_Off, max_disp_Off, RMS_On, med_disp_On, max_disp_On])
        
            header = 'RMS_Off med_disp_Off max_disp_Off RMS_On med_disp_On max_disp_On'
            if not os.path.exists(out_dir+'MotionMetrics_'+mot+'/'):
                os.makedirs(out_dir+'MotionMetrics_'+mot+'/')
            np.savetxt(out_dir+'MotionMetrics_'+mot+'/'+sequ+save+'.txt', save_list, header=header)






''' (2) load all data and calculate p-values '''

# stack p-values first STILL, then NOD, then SHAKE, in each cathegory: RMS, then median, then maximum
p_values = []
effect_size = []
for mot in ['STILL', 'NOD', 'SHAKE']:
    RMS, median, maxim, descr = [], [], [], []
    for sequ in sequs:
        if mot == 'SHAKE' and sequ != 'T1_MPR':
            continue
        # find the most recent data:
        file = glob.glob(out_dir+'MotionMetrics_'+mot+'/'+sequ+'*.txt')
        if len(file)>0:
            file = SortFiles(file)
            
        metrics = np.loadtxt(file[0], unpack=True, skiprows=1)
        RMS.append(metrics[0])
        RMS.append(metrics[3])
        median.append(metrics[1])
        median.append(metrics[4])
        maxim.append(metrics[2])
        maxim.append(metrics[5])
        descr.append(sequ)
        descr.append(sequ)
        
    p_val, ind, alt, ES = PerformWilcoxonMotion(['RMS', 'Med', 'Max'], mot, [RMS, median, maxim], out_dir, save)
    p_values.append(p_val)
    effect_size.append(ES)
        
#correct for multiple comparisons:
p_values = np.concatenate(p_values)
rej, p_values_cor, _, __ = multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

p_values_cor_still = p_values_cor[0:18]
p_values_cor_nod = p_values_cor[18:36]
p_values_cor_shake = p_values_cor[36:]
# check that everything correct!!!
ind_p = np.array([[0,1], [2,3], [4,5], [6,7], [8,9], [10,11]])
ind_sh = ind = np.array([[0,1]])



''' (3) Load the latest data and plot: '''

for mot in ['STILL', 'NOD']:
    RMS, median, maxim, descr = [], [], [], []
    for sequ in sequs:
        # find the most recent data:
        file = glob.glob(out_dir+'MotionMetrics_'+mot+'/'+sequ+'*.txt')
        if len(file)>0:
            file = SortFiles(file)
            
        metrics = np.loadtxt(file[0], unpack=True, skiprows=1)
        RMS.append(metrics[0])
        RMS.append(metrics[3])
        median.append(metrics[1])
        median.append(metrics[4])
        maxim.append(metrics[2])
        maxim.append(metrics[5])
        descr.append(sequ)
        descr.append(sequ)
    
    if mot == 'NOD':
        mot_s = 'SHAKE'
        RMS_s, median_s, maxim_s, descr_s = [], [], [], []
        sequ = 'T1_MPR'
        # find the most recent data:
        file = glob.glob(out_dir+'MotionMetrics_'+mot_s+'/'+sequ+'*.txt')
        if len(file)>0:
            file = SortFiles(file)
            
        metrics = np.loadtxt(file[0], unpack=True, skiprows=1)
        RMS_s.append(metrics[0])
        RMS_s.append(metrics[3])
        median_s.append(metrics[1])
        median_s.append(metrics[4])
        maxim_s.append(metrics[2])
        maxim_s.append(metrics[5])
        descr_s.append(sequ)
        descr_s.append(sequ)
        mean_RMS_s, mean_med_s, mean_max_s = [], [], []
        for i in range(len(RMS_s)):
            mean_RMS_s.append(np.mean(RMS_s[i]))
            mean_med_s.append(np.mean(median_s[i]))
            mean_max_s.append(np.mean(maxim_s[i]))
        
    
    mean_RMS, mean_med, mean_max = [], [], []
    for i in range(len(RMS)):
        mean_RMS.append(np.mean(RMS[i]))
        mean_med.append(np.mean(median[i]))
        mean_max.append(np.mean(maxim[i]))
    colors = []
    labels = []
    for i in range(int(len(RMS)/2)):
        colors.append('tab:orange')
        colors.append('tab:blue')
        labels.append('MOCO_OFF')
        labels.append('MOCO_ON')
    small = dict(markersize=3)
        
    
    # change 'TIRM' into STIR for descr:    
    # change 'TRACEW_B0' to 'DWI'
    for i in range(len(descr)):
        if descr[i]=='TRACEW_B0':
            descr[i]='DWI'
        if descr[i]=='T1_TIRM':
            descr[i]='T1_STIR'
        if descr[i]=='T2STAR':
            descr[i]='T2*'
    
    if mot=='STILL':
        plt.figure(figsize=(10,10))
        plt.subplot(3,1,1)
    else:
        plt.figure(figsize=(13, 11.5))

        ax1=plt.subplot2grid((3,5), (0,0), colspan=4)
    MakeBoxplot(RMS, colors)
    for i in range(len(mean_RMS)):
        plt.plot(i+1, mean_RMS[i], '.', c=colors[i], ls='')
    if mot=='STILL':
        for y1, y2 in zip(RMS[2], RMS[3]):
                plt.plot([3,4], [y1, y2], 'gray', lw=0.7)
        for y1, y2 in zip(RMS[8], RMS[9]):
                plt.plot([9,10], [y1, y2], 'gray', lw=0.7)
    plt.title(mot+' scans', fontsize=17)
    plt.ylabel('RMS displacement [mm]')
    plt.xticks(ticks=np.arange(1, len(RMS)+1), labels=descr)
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.2+lim[0])
    if mot=='STILL':
        use = 0
        low = -1
    if mot == 'NOD':
        use = 1
        low = -4
    
    max_RMS = [np.amax(RMS[i]) for i in range(len(RMS))]
    if mot=='STILL':
        Show_Stars(p_values_cor_still[0:6], ind_p, np.arange(1, len(RMS)+1), max_RMS)
        
        plt.subplot(3,1,2)
    else:
        Show_Stars(p_values_cor_nod[0:6], ind_p, np.arange(1, len(RMS)+1), max_RMS)
        
        ax2=plt.subplot2grid((3,5), (1,0), colspan=4)
    
    MakeBoxplot(median, colors)
    for i in range(len(mean_RMS)):
        plt.plot(i+1, mean_med[i], '.', c=colors[i], ls='')
    if mot=='STILL':
        for y1, y2 in zip(median[2], median[3]):
                plt.plot([3,4], [y1, y2], 'gray', lw=0.7)
        for y1, y2 in zip(median[8], median[9]):
                plt.plot([9,10], [y1, y2], 'gray', lw=0.7)
    plt.ylabel('Median displacement [mm]')
    plt.xticks(ticks=np.arange(1, len(RMS)+1), labels=descr)
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.2+lim[0])
    
    max_med = [np.amax(median[i]) for i in range(len(RMS))]
    if mot=='STILL':
        Show_Stars(p_values_cor_still[6:12], ind_p, np.arange(1, len(RMS)+1), max_med)
        plt.subplot(3,1,3)
    else:
        Show_Stars(p_values_cor_nod[6:12], ind_p, np.arange(1, len(RMS)+1), max_med)
        
        ax3=plt.subplot2grid((3,5), (2,0), colspan=4)
    MakeBoxplot(maxim, colors)
    for i in range(0,2):
        plt.plot(i+1, mean_max[i], '.', c=colors[i], ls='', label=labels[i])
    for i in range(2,len(mean_RMS)):
        plt.plot(i+1, mean_max[i], '.', c=colors[i], ls='')
    if mot=='STILL':
        for y1, y2 in zip(maxim[2], maxim[3]):
                plt.plot([3,4], [y1, y2], 'gray', lw=0.7)
        for y1, y2 in zip(maxim[8], maxim[9]):
                plt.plot([9,10], [y1, y2], 'gray', lw=0.7)
    plt.ylabel('Maximum displacement [mm]')
    plt.xticks(ticks=np.arange(1, len(RMS)+1), labels=descr)
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
    if mot=='STILL':
        use = 0
        low = -3
    
    max_max = [np.amax(maxim[i]) for i in range(len(RMS))]
    if mot=='STILL':
        Show_Stars(p_values_cor_still[12:], ind_p, np.arange(1, len(RMS)+1), max_max)
    else:
        Show_Stars(p_values_cor_nod[12:], ind_p, np.arange(1, len(RMS)+1), max_max)
    
    plt.legend(loc='upper center', ncol = 2, bbox_to_anchor=(0.5, -0.3), fontsize=11)
    
    if mot == 'NOD':
        ax4=plt.subplot2grid((3,5), (0,4))
        ax1.get_shared_y_axes().join(ax1, ax4)
        MakeBoxplot(RMS_s, colors)
        for i in range(len(mean_RMS_s)):
            plt.plot(i+1, mean_RMS_s[i], '.', c=colors[i], ls='')
        plt.title(mot_s+' scans', fontsize=17)
        plt.xticks(ticks=np.arange(1, len(RMS_s)+1), labels=descr_s)
        lim = plt.gca().get_ylim()
        plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
        max_RMS = [np.amax(RMS_s[i]) for i in range(len(RMS_s))]
        Show_Stars(np.array([p_values_cor_shake[0]]), ind_sh, np.arange(1, len(RMS)+1), max_RMS)
        
        ax5=plt.subplot2grid((3,5), (1,4))
        ax2.get_shared_y_axes().join(ax2, ax5)
        MakeBoxplot(median_s, colors)
        for i in range(len(mean_RMS_s)):
            plt.plot(i+1, mean_med_s[i], '.', c=colors[i], ls='')
        plt.xticks(ticks=np.arange(1, len(RMS_s)+1), labels=descr_s)
        lim = plt.gca().get_ylim()
        plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
        max_med = [np.amax(median_s[i]) for i in range(len(RMS_s))]
        Show_Stars(np.array([p_values_cor_shake[1]]), ind_sh, np.arange(1, len(RMS)+1), max_med)
        
        ax6=plt.subplot2grid((3,5), (2,4))
        ax3.get_shared_y_axes().join(ax3, ax6)
        MakeBoxplot(maxim_s, colors)
        for i in range(0,2):
            plt.plot(i+1, mean_max_s[i], '.', c=colors[i], ls='', label=labels[i])
        for i in range(2,len(mean_RMS_s)):
            plt.plot(i+1, mean_max_s[i], '.', c=colors[i], ls='')
        plt.xticks(ticks=np.arange(1, len(RMS_s)+1), labels=descr_s)
        lim = plt.gca().get_ylim()
        plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
        max_max = [np.amax(maxim[i]) for i in range(len(RMS_s))]
        Show_Stars(np.array([p_values_cor_shake[2]]), ind_sh, np.arange(1, len(RMS)+1), max_max)
           
    plt.subplots_adjust(hspace=0.3, wspace=0.4)
    plt.savefig(out_dir+'Boxplot_'+mot+save, bbox_inches='tight', dpi=200)
    plt.show()
    
    
    # print statistic over all sequences:
    print(mot)
    tmp = np.concatenate(RMS).ravel()
    print('RMS: ', np.median(tmp), np.std(tmp))
    tmp = np.concatenate(median).ravel()
    print('Median: ', np.median(tmp), np.std(tmp))
    tmp = np.concatenate(maxim).ravel()
    print('Maximum: ', np.median(tmp), np.std(tmp))
    
            
            
            
            
            
        
        
