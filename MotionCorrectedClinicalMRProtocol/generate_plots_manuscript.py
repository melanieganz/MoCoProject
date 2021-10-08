import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import string
import nibabel as nib
from statsmodels.stats.multitest import multipletests
from utils import Show_Stars, SortFiles, MakeBoxplot, DrawLines2, DrawLines
from statistical_tests import PerformWilcoxonMotion, PerformWilcoxonAllImg
plt.style.use('seaborn-whitegrid')
matplotlib.rc('axes', edgecolor='black')
from matplotlib.ticker import ScalarFormatter

out_dir = '/data1/hannah/Thesis_paper/Plots/'
in_dir = '/data1/hannah/Motion_Estimates/Comparison/'

subdir = [] 
for i in range(1,10):
    subdir.append('HC_0'+str(i)+'/')
for i in range(10,20):
    subdir.append('HC_'+str(i)+'/')
for i in range(20,23):
    subdir.append('HC_'+str(i)+'/')
sequs = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_TIRM', 'T2STAR', 'DIFF']

save = '_10_08'


''' (1) Plot motion data: '''
plot_motion = False

if plot_motion:
    
    # calculate p-values - first STILL, then NOD, then SHAKE, in each 
    # cathegory: RMS, then median, then maximum
    p_values = []
    effect_size = []
    for mot in ['STILL', 'NOD', 'SHAKE']:
        RMS, median, maxim, descr = [], [], [], []
        for sequ in sequs:
            if mot == 'SHAKE' and sequ != 'T1_MPR':
                continue
            # find the most recent data:
            file = glob.glob(in_dir+'MotionMetrics_'+mot+'/'+sequ+'*.txt')
            file = [f for f in file if 'mid' not in f]  #sort out 'mid'
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
            
        p_val, ind, alt, ES = PerformWilcoxonMotion(['RMS', 'Med', 'Max'], mot, 
                                                    [RMS, median, maxim], 
                                                    in_dir, save,)
        p_values.append(p_val)
        effect_size.append(ES)
            
    #correct for multiple comparisons:
    p_values = np.concatenate(p_values)
    rej, p_values_cor, _, __ = multipletests(p_values, alpha=0.05, 
                                             method='fdr_bh', is_sorted=False, 
                                             returnsorted=False)
    
    p_values_cor_still = p_values_cor[0:18]
    p_values_cor_nod = p_values_cor[18:36]
    p_values_cor_shake = p_values_cor[36:]
    ind_p = np.array([[0,1], [2,3], [4,5], [6,7], [8,9], [10,11]])
    ind_sh = np.array([[0,1]])
    
    
    # plot:
    num = 0
    plt.figure(figsize=(13,10))
    for mot in ['STILL', 'NOD']:
        RMS, median, maxim, descr = [], [], [], []
        for sequ in sequs:
            # find the most recent data:
            file = glob.glob(in_dir+'MotionMetrics_'+mot+'/'+sequ+'*.txt')
            file = [f for f in file if 'mid' not in f]  # sort out mid
            if len(file)>0:
                file = SortFiles(file)
                
            # load the most recent file:
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
            file = glob.glob(in_dir+'MotionMetrics_'+mot_s+'/'+sequ+'*.txt')
            file = [f for f in file if 'mid' not in f] # sort out mid
            if len(file)>0:
                file = SortFiles(file)
                
            # load the most recent file
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
          
        # calculate mean values and define colors
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
            labels.append('without PMC')
            labels.append('with PMC')
        small = dict(markersize=3)
            
        # change 'TIRM' into STIR for descr:    
        # and change 'TRACEW_B0' to 'DWI'
        for i in range(len(descr)):
            if descr[i]=='TRACEW_B0':
                descr[i]='DWI'
            if descr[i]=='T1_TIRM':
                descr[i]='T1_STIR'
            if descr[i]=='T2STAR':
                descr[i]='T2*'
        
        # plot the data:
        if mot=='STILL':
            ax0=plt.subplot2grid((4,5), (0,0), colspan=4)
            ax = ax0
        else:
            ax1=plt.subplot2grid((4,5), (2,0), colspan=4)
            ax = ax1
        MakeBoxplot(RMS, colors)
        for i in range(len(mean_RMS)):
            plt.plot(i+1, mean_RMS[i], '.', c=colors[i], ls='')
        if mot=='STILL':
            for y1, y2 in zip(RMS[4], RMS[5]):
                    plt.plot([5,6], [y1, y2], 'darkslategray', lw=0.7)
            for y1, y2 in zip(RMS[8], RMS[9]):
                    plt.plot([9,10], [y1, y2], 'darkslategray', lw=0.7)

        plt.title(mot+' scans', fontsize=17)
        plt.ylabel('RMS displ [mm]')
        plt.xticks(ticks=np.arange(1, len(RMS)+1), labels=descr)
        lim = plt.gca().get_ylim()
        plt.ylim(lim[0],(lim[1]-lim[0])*1.2+lim[0])
        ax.text(-0.1, .95, string.ascii_lowercase[num], transform=ax.transAxes,
                size=20, weight='bold')
        num += 1
        if mot=='STILL':
            use = 0
            low = -1
        if mot == 'NOD':
            use = 1
            low = -4
        
        max_RMS = [np.amax(RMS[i]) for i in range(len(RMS))]
        if mot=='STILL':
            Show_Stars(p_values_cor_still[0:6], ind_p, np.arange(1, len(RMS)+1), 
                       max_RMS, dh=0.05, col='black')
            
            ax01=plt.subplot2grid((4,5), (1,0), colspan=4)
            ax = ax01
        else:
            Show_Stars(p_values_cor_nod[0:6], ind_p, np.arange(1, len(RMS)+1), 
                       max_RMS, col='black')
            
            ax2=plt.subplot2grid((4,5), (3,0), colspan=4)
            ax = ax2
        MakeBoxplot(maxim, colors)
        for i in range(0,2):
            plt.plot(i+1, mean_max[i], '.', c=colors[i], ls='', 
                     label=labels[i])
        for i in range(2,len(mean_RMS)):
            plt.plot(i+1, mean_max[i], '.', c=colors[i], ls='')
        if mot=='STILL':
            for y1, y2 in zip(maxim[4], maxim[5]):
                    plt.plot([5,6], [y1, y2], 'darkslategray', lw=0.7)
            for y1, y2 in zip(maxim[8], maxim[9]):
                    plt.plot([9,10], [y1, y2], 'darkslategray', lw=0.7)
        plt.ylabel('Maximum displ [mm]')
        plt.xticks(ticks=np.arange(1, len(RMS)+1), labels=descr)
        lim = plt.gca().get_ylim()
        plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
        ax.text(-0.1, .95, string.ascii_lowercase[num], transform=ax.transAxes,
                size=20, weight='bold')
        num += 1
        if mot=='STILL':
            use = 0
            low = -3       
        max_max = [np.amax(maxim[i]) for i in range(len(RMS))]
        if mot=='STILL':
            Show_Stars(p_values_cor_still[12:], ind_p, np.arange(1, len(RMS)+1), 
                       max_max, col='black')
            legend = plt.legend(loc='upper center', bbox_to_anchor=(1.16, 1.3), 
                                fontsize=13, frameon=True)
        else:
            Show_Stars(p_values_cor_nod[12:], ind_p, np.arange(1, len(RMS)+1), 
                       max_max, col='black')
        
        if mot == 'NOD':
            ax4=plt.subplot2grid((4,5), (2,4))
            ax = ax4
            ax1.get_shared_y_axes().join(ax1, ax4)
            MakeBoxplot(RMS_s, colors)
            for i in range(len(mean_RMS_s)):
                plt.plot(i+1, mean_RMS_s[i], '.', c=colors[i], ls='')
            plt.title(mot_s+' scans', fontsize=17)
            plt.xticks(ticks=np.arange(1, len(RMS_s)+1), labels=descr_s)
            lim = plt.gca().get_ylim()
            plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
            ax.text(-0.3, .95, string.ascii_lowercase[num], 
                    transform=ax.transAxes, size=20, weight='bold')
            num += 1
            max_RMS = [np.amax(RMS_s[i]) for i in range(len(RMS_s))]
            Show_Stars(np.array([p_values_cor_shake[0]]), ind_sh, 
                       np.arange(1, len(RMS)+1), max_RMS, col='black')
            
            ax5=plt.subplot2grid((4,5), (3,4))
            ax = ax5
            ax2.get_shared_y_axes().join(ax2, ax5)
            
            MakeBoxplot(maxim_s, colors)
            for i in range(0,2):
                plt.plot(i+1, mean_max_s[i], '.', c=colors[i], ls='', 
                         label=labels[i])
            for i in range(2,len(mean_RMS_s)):
                plt.plot(i+1, mean_max_s[i], '.', c=colors[i], ls='')
            plt.xticks(ticks=np.arange(1, len(RMS_s)+1), labels=descr_s)
            lim = plt.gca().get_ylim()
            plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
            ax.text(-0.3, .95, string.ascii_lowercase[num], 
                    transform=ax.transAxes, size=20, weight='bold')
            num += 1
            max_max = [np.amax(maxim[i]) for i in range(len(RMS_s))]
            Show_Stars(np.array([p_values_cor_shake[2]]), ind_sh, 
                       np.arange(1, len(RMS)+1), max_max, col='black')
    
    legend.get_frame().set_linewidth(2)           
    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    plt.savefig(out_dir+'Boxplot_motion_'+save+'.tiff', format='tiff', 
                bbox_inches='tight', dpi=200)
    plt.savefig(out_dir+'Boxplot_motion_'+save+'.png', format='png', 
                bbox_inches='tight', dpi=200)
    plt.show()



''' (2) Import image quality metrics: '''
in_dir = '/data1/hannah/Metrics_Results/'
out_dir_metrics = '/data1/hannah/Metrics_Results/Comparison/'
withRR = True
onlyRR = False    
quality_scores = True
show_stat_test = True

sequs = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_TIRM', 'T2STAR', 'TRACEW_B1000']

ssims, psnrs, tgs, qss, names_tes = [], [], [], [], []
p_ssims, p_psnrs, p_tgs, p_qss, ind_ps = [], [], [], [], []
for sequ in sequs:
    n = len(subdir)
    if sequ == 'T1_MPR':
        m = 10
        names = np.array(['MOCO_ON_STILL', 'MOCO_ON_SHAKE_RR', 'MOCO_ON_SHAKE', 
                          'MOCO_ON_NOD_RR', 'MOCO_ON_NOD', 'MOCO_OFF_STILL', 
                          'MOCO_OFF_SHAKE_RR', 'MOCO_OFF_SHAKE', 
                          'MOCO_OFF_NOD_RR', 'MOCO_OFF_NOD'])
    else:
        m = 6
        names = np.array(['MOCO_ON_STILL', 'MOCO_ON_NOD_RR', 'MOCO_ON_NOD', 
                          'MOCO_OFF_STILL', 'MOCO_OFF_NOD_RR', 'MOCO_OFF_NOD'])
    ssim, psnr, tg_, qs = np.zeros((n,m)), np.zeros((n,m)), np.zeros((n,m)),  np.zeros((n,m))
    i=0
    print('Files used for calculation:')
    for sub in subdir:
        # for DWI, T2FLAIR and T2STAR not all volunteers available:
        if sequ in ['ADC', 'TRACEW_B0', 'TRACEW_B1000', 'T2_FLAIR'] and sub in ['HC_0'+str(i)+'/' for i in range(1,10)]:
            continue
        if sequ in ['ADC', 'TRACEW_B0', 'TRACEW_B1000', 'T2_FLAIR'] and sub in ['HC_'+str(i)+'/' for i in range(10,13)]:
            continue
        if sequ == 'T2STAR' and sub in ['HC_20/', 'HC_21/', 'HC_22/']:
            continue
        
        folder = in_dir+sub  
        file_ = glob.glob(folder+'Values_*'+sequ+'*')
        file = [f for f in file_ if f.find('Retro')==-1]
        if len(file)>0:
            # get the most recent file:
            file = SortFiles(file)
            print(file[0])
            tmp1 = np.loadtxt(file[0], unpack=True, dtype=str, usecols=0, 
                              skiprows=1)
            tmp2, tmp3, tmp4 = np.loadtxt(file[0], unpack=True, 
                                          usecols=(1,2,3), skiprows=1)
            incl = []
            for na,s,p,t in zip(tmp1, tmp2, tmp3, tmp4):
                if sequ in ['T2STAR', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
                    if 'NOD' in na:
                        na = na+'_RR'
                
                ind = np.where(names==na)[0][0]
                ssim[i,ind], psnr[i,ind], tg_[i,ind] = s, p, t
                incl.append(ind)
                # 0 entires represent non-existing scans
            
            # do the same for quality scores:
            if sequ not in ['T2STAR', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
                file = glob.glob(folder+'QualityScores*'+sequ+'*')
                if len(file)>0:
                    # get the most recent file:
                    file = SortFiles(file)
                    print(file[0])
                    tmp1 = np.loadtxt(file[0], unpack=True, dtype=str, 
                                      usecols=0, skiprows=1)
                    tmp2 = np.loadtxt(file[0], unpack=True, usecols=1, 
                                      skiprows=1)
                    incl = []
                    for na,q in zip(tmp1, tmp2):
                        if 'RETRO' not in na:
                            ind = np.where(names==na)[0][0]
                            qs[i,ind] = q
                            incl.append(ind)
                            
            # for DWI quality ranks instead of quality scores:                
            elif sequ == 'TRACEW_B1000':
                file = glob.glob(folder+'Ranks_DWI**') 
                if len(file)>0:
                    # get the most recent file:
                    file = SortFiles(file)
                    print(file[0])
                    tmp1 = np.loadtxt(file[0], unpack=True, dtype=str, 
                                      usecols=0, skiprows=1)
                    tmp2 = np.loadtxt(file[0], unpack=True, usecols=1, 
                                      skiprows=1)
                    incl = []
                    for na,q in zip(tmp1, tmp2):
                        if 'NOD' in na:
                            na = na+'_RR'
                        if 'RETRO' not in na:
                            ind = np.where(names==na)[0][0]
                            qs[i,ind] = q
                            incl.append(ind)
            else:
                qs = np.zeros_like(ssim)
            i+=1

    names = np.tile(names, (n,1))

    # check that names are the same for all volunteers:
    for i in range(1, len(subdir)):
        if np.all(names[0]==names[i], axis=0)==False:
            print('ERROR: Names of 0 and '+str(i)+' are not the same!')

    names = names[0]

    ''' sort out values == 0 (corresponding to missing scans) '''
    ssim[ssim==0] = np.nan
    psnr[psnr==0] = np.nan
    tg_[tg_==0] = np.nan
    qs[qs==0] = np.nan


    ''' perform parametric tests:'''
    # resort names and data (only for test), final resorting will be performed
    # after RR scans are potentially sorted out:
    names_ch = []
    for n in names:
        ind=n.find('_', 5)
        if 'OFF' in n:
            tmp = 'C'
        elif 'ON' in n:
            tmp = 'A'
        elif 'RETRO' in n:
            tmp = 'B'
        names_ch.append(n[ind+1:]+'_'+tmp)

    names_ch = np.array(names_ch)
    ind = np.argsort(names_ch)[::-1]

    ssim_p, psnr_p, tg_p, qs_p = ssim[:,ind], psnr[:,ind], tg_[:,ind], qs[:,ind]
    names_p = names_ch[ind]

    p_ssim, rej_ssim, ind_p, alt = PerformWilcoxonAllImg('SSIM', ssim_p, sequ, 
                                                         out_dir_metrics, save)
    p_psnr, rej_psnr, ind_p, alt = PerformWilcoxonAllImg('PSNR', psnr_p, sequ, 
                                                         out_dir_metrics, save)
    p_tg, rej_tg, ind_p, alt = PerformWilcoxonAllImg('TG', tg_p, sequ, 
                                                     out_dir_metrics, save)
    p_qs, rej_qs, ind_p, alt = PerformWilcoxonAllImg('QS', qs_p, sequ, 
                                                     out_dir_metrics, save)

    # sort out the relevant tests:
    if onlyRR == True:
        if sequ == 'T1_MPR':
            rel_tests = [0,1,4]
            p_ssim, p_psnr = p_ssim[rel_tests], p_psnr[rel_tests]
            p_tg, p_qs = p_tg[rel_tests], p_qs[rel_tests]
            ind_p = np.array([[0,1], [2,3], [4,5]])
        else:
            rel_tests = [0,1]
            p_ssim, p_psnr = p_ssim[rel_tests], p_psnr[rel_tests]
            p_tg, p_qs = p_tg[rel_tests], p_qs[rel_tests]
            ind_p = np.array([[0,1], [2,3]])

    if withRR == False:
            findRRs = np.char.find(names, 'RR')
            findRRs = (findRRs<0)
            names = names[findRRs]
            ssim, psnr, tg_ = ssim[findRRs], psnr[findRRs], tg_[findRRs]
            ssim = np.reshape(ssim, (len(subdir), int(len(ssim)/len(subdir))))
            psnr = np.reshape(psnr, (len(subdir), int(len(psnr)/len(subdir))))
            tg_ = np.reshape(tg_, (len(subdir), int(len(tg_)/len(subdir))))
            qs = np.reshape(qs, (len(subdir), int(len(qs)/len(subdir))))
            names = np.reshape(names, (len(subdir), int(len(names)/len(subdir))))

    if onlyRR == True:
        findRRs = np.char.find(names, 'RR')
        findstills = np.char.find(names, 'STILL')
        both = findRRs*findstills
        both = (both<0)
        names = names[both]
        ssim, psnr, tg_, qs = ssim[:,both], psnr[:,both], tg_[:,both], qs[:,both]


    std_ssim, std_psnr  = np.nanstd(ssim, axis=0), np.nanstd(psnr, axis=0) 
    std_tg_, std_qs = np.nanstd(tg_, axis=0), np.nanstd(qs, axis=0)
    mean_ssim, mean_psnr = np.nanmean(ssim, axis=0), np.nanmean(psnr, axis=0) 
    mean_tg_, mean_qs = np.nanmean(tg_, axis=0), np.nanmean(qs, axis=0)

    # Final resorting:
    names_ch = []
    for n in names:
        ind=n.find('_', 5)
        if 'OFF' in n:
            tmp = 'C'
        elif 'ON' in n:
            tmp = 'A'
        elif 'RETRO' in n:
            tmp = 'B'
        names_ch.append(n[ind+1:]+'_'+tmp)

    names_ch = np.array(names_ch)
    ind = np.argsort(names_ch)[::-1]
    ssim, psnr, tg_, qs = ssim[:,ind], psnr[:,ind], tg_[:,ind], qs[:,ind]
    names = names[ind]
    mean_ssim, mean_psnr =  mean_ssim[ind], mean_psnr[ind]
    mean_tg_, mean_qs = mean_tg_[ind], mean_qs[ind]

    names_te = []
    for i in range(int(len(names))):
        tmp = names[i]
        if 'RETRO' in tmp:
            app = 'retrospective MoCo'
        elif 'ON' in tmp:
            app = 'prospective MoCo'
        elif 'OFF' in tmp:
            app = 'no MoCo'
        if onlyRR == False:
            if 'RR' in tmp or 'STILL' in tmp:
                app = app + ' no REAC'
            else:
                app = app + ' with REAC'
        names_te.append(app)

    
    # sort out nans (whole rows for subjects wihtout that sequence)
    mask = np.where(np.isnan(ssim[:,0])==False)[0]
    ssim = ssim[mask]
    mask = np.where(np.isnan(psnr[:,0])==False)[0]
    psnr = psnr[mask]
    mask = np.where(np.isnan(tg_[:,0])==False)[0]
    tg_ = tg_[mask]
    mask = np.where(np.isnan(qs[:,0])==False)[0]
    qs = qs[mask]
    
    # normalise tg values by dividing with mean value of off still
    val = np.mean(tg_[:,0])
    tg_ = tg_/val
    

    # new color code:
    color_dict = {'retrospective MoCo': 'tab:green', 
                  'prospective MoCo': 'tab:blue', 'no MoCo':'tab:orange'}
    if onlyRR == False:
        color_dict = {'retrospective MoCo no REAC': 'tab:green', 
                      'prospective MoCo no REAC': 'tab:blue', 
                      'no MoCo no REAC':'tab:orange',
                      'retrospective MoCo with REAC': 'tab:gray', 
                      'prospective MoCo with REAC': 'tab:olive', 
                      'no MoCo with REAC':'tab:cyan'}

    ssims.append(ssim)
    psnrs.append(psnr)
    tgs.append(tg_)
    qss.append(qs)
    names_tes.append(names_te)
    p_ssims.append(p_ssim)
    p_psnrs.append(p_psnr)
    p_tgs.append(p_tg)
    p_qss.append(p_qs)
    ind_ps.append(ind_p)



''' (3) Plot image quality metrics for still scans: '''
plot_still = True
if plot_still:
    metrics = [tgs, qss, np.array([qss[-1]])]
    p_values = [p_tgs, p_qss[0:4], np.array([p_qss[-1]])]
    labels = ['Tenengrad', 'Observer Scores', 'Quality Rank']
    colors = ['tab:orange', 'tab:blue', 'tab:orange', 'tab:blue', 'tab:orange', 
              'tab:blue', 'tab:orange', 'tab:blue', 'tab:orange', 'tab:blue', 
              'tab:orange', 'tab:blue']
    fig_labels = ['without PMC', 'with PMC']
    small = dict(markersize=3)
    x = np.arange(1,13)
    a = [0,2,4,6,8,10]
    b = [1,3,5,7,9,11]
    c = [1,3,5,7,9,11]
    d = [2,4,6,8,10,12]
    
    num = 0
    plt.figure(figsize=(10,9))
    for i, metric, p, lab in zip(range(0,3), metrics, p_values, labels):
        m_still = []
        means = []
        for m in metric:
            if len(m)>0:
                m_still.append(m[:,0])
                m_still.append(m[:,1])
                means.append(np.nanmean(m[:,0]))
                means.append(np.nanmean(m[:,1]))
        p_still = []
        for pval in p:
            p_still.append(pval[0])
        
        if i == 0:
            ax = plt.subplot2grid((2,6), (0,0), colspan=6)
            box1 = plt.boxplot(m_still, flierprops=small)
            for j in range(len(means)):
                plt.errorbar(x[j], means[j], yerr=None, color=colors[j], fmt='.', capsize=3)
            ticklabels = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_STIR', 'T2*', 'TRACEW']
            ticks = [1.5, 3.5, 5.5, 7.5, 9.5, 11.5]
            plt.xticks(labels=ticklabels, ticks=ticks, fontsize=14)
                
        elif i == 1:
            ax = plt.subplot2grid((2,6), (1,0), colspan=4)
            box1 = plt.boxplot(m_still[0:8], flierprops=small)
            for j in range(len(means)-2):
                plt.errorbar(x[j], means[j], yerr=None, color=colors[j], fmt='.', capsize=3)
            for j in range(0,2):
                plt.errorbar(x[j], means[j], yerr=None, color=colors[j], fmt='.', 
                             capsize=3, label=fig_labels[j])
            ticklabels = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_STIR']
            ticks = [1.5, 3.5, 5.5, 7.5]
            plt.xticks(labels=ticklabels, ticks=ticks, fontsize=14)
            
        else:
            ax = plt.subplot2grid((2,6), (1,5), colspan=1)
            box1 = plt.boxplot(m_still, flierprops=small)
            for j in range(len(means)):
                plt.errorbar(x[j], means[j], yerr=None, color=colors[j], fmt='.', capsize=3)
            ticklabels = ['DWI']
            ticks = [1.5]
            plt.xticks(labels=ticklabels, ticks=ticks, fontsize=14)
        
        for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
            patch.set(color=color, lw=1.7)
            patch2.set(color='k', lw=1.7)
           
        plt.ylabel(lab, fontsize=15)
        if show_stat_test == True:
            maxi = []
            for v in m_still:
                maxi.append(np.amax(v))
            if i == 0:
                indices = [[0,1], [2,3], [4,5], [6,7], [8,9], [10,11]]
                Show_Stars(np.array(p_still), indices[0:len(p_still)], x, maxi, 
                           col='black')
            elif i ==1:
                indices = [[0,1], [2,3], [4,5], [6,7]]
                Show_Stars(np.array(p_still), indices[0:len(p_still)], x[0:8], maxi[0:8], 
                           col='black')
            else:
                indices = [[0,1]]
                Show_Stars(np.array(p_still), indices[0:len(p_still)], x, maxi, 
                           col='black')
            lim = plt.gca().get_ylim()
            plt.ylim(lim[0],(lim[1]-lim[0])*1.05+lim[0])
        DrawLines2(a[0:len(p_still)],b[0:len(p_still)],c[0:len(p_still)],
                   d[0:len(p_still)],m_still, lw=1, col='darkslategray')
        
        if i == 2:
            ax.text(-0.6, 0.95, string.ascii_lowercase[num], 
                    transform=ax.transAxes, size=24, weight='bold')
        else:
            ax.text(-0.1, 0.95, string.ascii_lowercase[num], 
                    transform=ax.transAxes, size=24, weight='bold')
        num += 1
       
        plt.yticks(fontsize=13)
        plt.tick_params('both', length=0)
        plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        if i == 0:
            xlim = plt.gca().get_xlim()
        if i == 1:
            legend =  plt.legend( loc='lower left', ncol=2, 
                                bbox_to_anchor=(0.4, -0.3), fontsize=14, 
                                frameon=True)
            plt.yticks(ticks=[2.5, 3, 3.5, 4, 4.5, 5])
        
    legend.get_frame().set_linewidth(2)           
    plt.subplots_adjust(hspace=0.2, wspace=0.3)
    plt.savefig(out_dir+'Metrics_still'+save+'.tiff', format='tiff', 
                bbox_inches='tight', dpi=200)
    plt.savefig(out_dir+'Metrics_still'+save+'.png', format='png', 
                bbox_inches='tight', dpi=200)
    plt.show()



''' (4) Plot image quality metrics for nodding scans: '''
plot_nod = True
if plot_nod:
    # first MPRAGE and FLAIR
    metrics = [ssims, psnrs, tgs, qss]
    p_values = [p_ssims, p_psnrs, p_tgs, p_qss[0:4]]
    labels = ['SSIM', 'PSNR', 'Tenengrad', 'Observer Scores']
    color_dict = {'prospective MoCo no REAC': 'tab:blue', 
                  'no MoCo no REAC':'tab:orange',
                  'prospective MoCo with REAC': 'tab:green', 
                  'no MoCo with REAC':'tab:cyan'}
    name_dict = {'prospective MoCo no REAC': 'with PMC without REAC',
                  'no MoCo no REAC':'without PMC without REAC',
                  'prospective MoCo with REAC': 'with PMC with REAC',
                  'no MoCo with REAC':'without PMC with REAC'}
    
    small = dict(markersize=3)
    x = np.concatenate((np.arange(1,9), np.arange(10,14)), axis=None)
    a = [0,2,4,6] 
    b = [1,3,5,7]
    c = [1,3,5,7]
    d = [2,4,6,8]
    a_ = [0,2]
    b_ = [1,3]
    c_ = [10,12]
    d_ = [11,13]
    
    num = 0
    plt.figure(figsize=(11,8))
    for i, metric, p, lab in zip(range(0,4), metrics, p_values, labels):
        means = []
        names = []
        # first MPR
        tmp = metric[0]
        m_mpr = tmp[:,2:]
        names = names_tes[0][2:]
        
        # then FLAIR
        tmp = metric[1]
        m_fl = tmp[:,2:]
        names = np.concatenate((names, names_tes[1][2:]), axis=None)
        
        p_mot = []
        indices = []
        for pval, ind in zip(p[0:2], ind_ps[0:2]):
            p_mot.append(pval[1:])
            indices.append(ind[1:])
            
        means = np.concatenate((np.nanmean(m_mpr, axis=0), 
                                np.nanmean(m_fl, axis=0)), axis=None)
        
        colors = []
        leg_labels = []
        for n in names:
            colors.append(color_dict[n])
            leg_labels.append(name_dict[n])
            
        ax = plt.subplot(2,2,i+1)
        box1 = plt.boxplot(m_mpr, flierprops=small, widths=0.5)
        for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
            patch.set(color=color, lw=1.5)
            patch2.set(color='k', lw=1.5)
        box2 = plt.boxplot(m_fl, positions=range(10,14), flierprops=small, 
                           widths=0.5)
        for patch, patch2, color in zip(box2['boxes'], box2['medians'], colors):
            patch.set(color=color, lw=1.5)
            patch2.set(color='k', lw=1.5)
        for j in range(len(means)):
            plt.errorbar(x[j], means[j], yerr=None, color=colors[j], fmt='.', 
                         capsize=3)
        if i == 2:
            for j in range(0,4):
                plt.errorbar(x[j], means[j], yerr=None, color=colors[j], 
                             fmt='.', capsize=3, label=leg_labels[j])
            legend = plt.legend( loc='lower left', ncol=2, 
                                bbox_to_anchor=(0.25, -0.4), fontsize=12, 
                                frameon=True)
        plt.ylabel(lab, fontsize=15)
        if show_stat_test == True:
            maxi = []
            for v in m_mpr.T:
                maxi.append(np.amax(v))
            Show_Stars(np.array(p_mot[0]), indices[0]-2, x, maxi, 
                       arange_dh='PAPER', col='black')
            for v in m_fl.T:
                maxi.append(np.amax(v))
            Show_Stars(np.array(p_mot[1]), indices[1]-2, x[8:], maxi, 
                       arange_dh='PAPER', col='black')
            lim = plt.gca().get_ylim()
            plt.ylim(lim[0],(lim[1]-lim[0])*1.01+lim[0])
        DrawLines(a,b,c,d,m_mpr, col='darkslategray')
        DrawLines(a_,b_,c_,d_,m_fl, col='darkslategray')
        ticklabels = ['T1_MPR\nSHAKE', 'T1_MPR\nNOD', 'T2_FLAIR\nNOD']
        ticks = [2.5, 6.5, 11.5]
        plt.xticks(labels=ticklabels, ticks=ticks, fontsize=12)
        ax.text(-0.18, 0.9, string.ascii_lowercase[num], transform=ax.transAxes,
                size=21, weight='bold')
        num += 1
        plt.yticks(fontsize=13)
        plt.tick_params('both', length=0)
        plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        if i == 3:
            plt.yticks(ticks=[1, 2, 3, 4, 5])
    
    legend.get_frame().set_linewidth(2)           
    plt.subplots_adjust(hspace=0.2, wspace=0.3)
    plt.savefig(out_dir+'Metrics_3D'+save+'.tiff', format='tiff', 
                bbox_inches='tight', dpi=200)
    plt.savefig(out_dir+'Metrics_3D'+save+'.png', format='png', 
                bbox_inches='tight', dpi=200)
    plt.show()
    
    
    
    # now TSE, STIR, T2* and TRACEW
    metrics = [ssims, psnrs, tgs, qss, np.array(qss[-1])]
    p_values = [p_ssims, p_psnrs, p_tgs, p_qss[0:4], np.array([p_qss[-1]])]
    labels = ['SSIM', 'PSNR', 'Tenengrad', 'Observer Scores', 'Quality Rank']
    color_dict = {'prospective MoCo no REAC': 'tab:blue', 
                  'no MoCo no REAC':'tab:orange',
                  'prospective MoCo with REAC': 'tab:green', 
                  'no MoCo with REAC':'tab:cyan'}
    name_dict = {'prospective MoCo no REAC': 'with PMC without REAC',
                  'no MoCo no REAC':'without PMC without REAC',
                  'prospective MoCo with REAC': 'with PMC with REAC',
                  'no MoCo with REAC':'without PMC with REAC'}   
    small = dict(markersize=3)
    
    num = 0
    plt.figure(figsize=(13,8))
    for i, metric, p, lab in zip(range(0,5), metrics, p_values, labels):
        means = []
        names = []
        m_nod = []
        
        if i == 4:
            tmp = metric
            m_nod= [tmp[:,2:4]]
            means = [np.nanmean(tmp[:,2:4], axis=0)]
            names = [names_tes[j][2:4]]
            
            p_mot = [p[0,1]]
            indices = [ind_ps[-1][1]]

        elif i in [0,1,2,3]:
            for j in range(2, 6):
                end = 6
                if j in [4,5]:
                    end = 4
                tmp = metric[j]
                m_nod.append(tmp[:,2:end])
                means.append(np.nanmean(tmp[:,2:end], axis=0))
                names.append(names_tes[j][2:end])
        
            p_mot = []
            indices = []
            for pval, ind in zip(p[2:6], ind_ps[2:6]):
                p_mot.append(pval[1:])
                indices.append(ind[1:])

        if i in [0,1,2]:
            ax = plt.subplot2grid((2,4), (i//2,i%2*2), colspan=2)
        if i == 3:
            ax = plt.subplot2grid((2,6), (1,3), colspan=2)
        if i == 4:
            ax = plt.subplot2grid((2,6), (1,5), colspan=1)
            
        N = 1
        lims = []
        
        for m, mean, p, index, name  in zip(m_nod, means, p_mot, indices, names):
            colors = []
            leg_labels = []
            for n in name:
                colors.append(color_dict[n])
                leg_labels.append(name_dict[n])
            if i in [0,1,2]:
                positions = np.arange(N, N+len(m.T))
                box1 = plt.boxplot(m, positions=positions, flierprops=small, 
                                   widths=0.5)
                for j in range(len(mean)):
                    plt.errorbar(positions[j], mean[j], yerr=None, color=colors[j], 
                                 fmt='.', capsize=3)          
                if i == 2 and N == 1:
                    for j in range(0,4):
                        plt.errorbar(positions[j], mean[j], yerr=None, 
                                     color=colors[j], fmt='.', 
                                     capsize=3, label=leg_labels[j])
                    legend = plt.legend( loc='lower left', ncol=2, 
                                        bbox_to_anchor=(0.45, -0.4), fontsize=12, 
                                        frameon=True)
                ticklabels = ['T2_TSE', 'T1_STIR', 'T2*', 'TRACEW']
                ticks = [2.5, 7.5, 11.5, 14.5]
                plt.xticks(labels=ticklabels, ticks=ticks, fontsize=12)
                    
            if i == 3:
                positions = np.arange(N, N+len(m.T))
                box1 = plt.boxplot(m[0:8], positions=positions[0:8], flierprops=small, 
                                   widths=0.5)
                for j in range(len(mean)):
                    plt.errorbar(positions[j], mean[j], yerr=None, color=colors[j], 
                                 fmt='.', capsize=3)
                ticklabels = ['T2_TSE', 'T1_STIR']
                ticks = [2.5, 7.5]
                plt.xticks(labels=ticklabels, ticks=ticks, fontsize=12)
                
            if i == 4:
                positions = np.arange(N, N+len(m.T))
                box1 = plt.boxplot(m, positions=positions, flierprops=small, 
                                   widths=0.5)
                for j in range(len(mean)):
                    plt.errorbar(positions[j], mean[j], yerr=None, color=colors[j], 
                                 fmt='.', capsize=3)
                ticklabels = ['DWI']
                ticks = [1.5]
                plt.xticks(labels=ticklabels, ticks=ticks, fontsize=12)
               
                
            for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
                patch.set(color=color, lw=1.5)
                patch2.set(color='k', lw=1.5)
            
            
            if show_stat_test == True:
                maxi = []
                for v in m.T:
                    maxi.append(np.amax(v))
                if i in [0,1,2]:
                    Show_Stars(np.array(p), index-2, positions, maxi, 
                               flexible_dh=True, col='black')
                if i == 3:
                    Show_Stars(np.array(p), index-2, positions[0:8], maxi[0:8], 
                               flexible_dh=True, col='black')
                
                if i == 4:
                    Show_Stars(np.array([p]), [index-2], positions, maxi, 
                               flexible_dh=True, col='black')
                
                lims.append(plt.gca().get_ylim())

            a = [0,2]
            b = [1,3]
            if N > 10 or i == 4:
                a = [0]
                b = [1]
                
            c = [positions[p] for p in a]
            d = [positions[p] for p in b]
            DrawLines(a,b,c,d,m, col='darkslategray')
            N = N+len(m.T)+1
            
        if show_stat_test:
            lims = np.array(lims)
            lim_min = np.amin(lims[:,0])
            lim_max = np.amax(lims[:,1])
            plt.ylim(lim_min,(lim_max-lim_min)*1.01+lim_min)
            
        plt.ylabel(lab, fontsize=15)      
        plt.yticks(fontsize=13)
        plt.tick_params('both', length=0)
        plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        if i == 4:
            ax.text(-0.51, 0.9, string.ascii_lowercase[num], transform=ax.transAxes,
                size=21, weight='bold')
        else:
            ax.text(-0.18, 0.9, string.ascii_lowercase[num], transform=ax.transAxes,
                size=21, weight='bold')                
        num += 1
        if i == 3:
            plt.yticks(ticks=[1, 2, 3, 4, 5])
        if i == 4:
            plt.yticks(ticks=[1, 2, 3, 4])
    
    legend.get_frame().set_linewidth(2)           
    plt.subplots_adjust(hspace=0.2, wspace=0.6)
    plt.savefig(out_dir+'Metrics_2D'+save+'.tiff', format='tiff', 
                bbox_inches='tight', dpi=200)
    plt.savefig(out_dir+'Metrics_2D'+save+'.png', format='png', 
                bbox_inches='tight', dpi=200)
    plt.show()



''' (5) Plot DWI analysis: '''
plot_DWI = False
if plot_DWI:
    # remember to add code for calculating the histograms somewhere!
    All_mean, All_std, All_sum = np.load(out_dir_metrics+'ADC_all_values_04_26.npy')
    All_mean = np.array([All_mean[:,0], All_mean[:,2], All_mean[:,1]]).T
    All_std = np.array([All_std[:,0], All_std[:,2], All_std[:,1]]).T
    All_sum = np.array([All_sum[:,0], All_sum[:,2], All_sum[:,1]]).T
    names_diff = ['prospective MoCo', 'no MoCo', 'prospective MoCo']
    
    p_mean, rej_mean, ind, altern = PerformWilcoxonAllImg('Mean', All_mean, 'ADC', 
                                                          out_dir_metrics, 
                                                          '_04_26', option='diff')
    p_std, rej_std, ind, altern = PerformWilcoxonAllImg('Std', All_std, 'ADC',
                                                        out_dir_metrics, 
                                                        '_04_26', option='diff')
    p_sum, rej_sum, ind, altern = PerformWilcoxonAllImg('Sum', All_sum, 'ADC', 
                                                        out_dir_metrics, 
                                                        '_04_26', option='diff')
    
    mean_mean = np.mean(All_mean, axis=0)
    mean_std = np.mean(All_std, axis=0)
    mean_sum = np.mean(All_sum, axis=0)
    
    color_dict = { 'prospective MoCo': 'tab:blue', 'no MoCo': 'tab:orange'}
    names = ['with PMC', 'without PMC', 'with PMC']
    
    colors_te = []
    for n in names_diff:
        colors_te.append(color_dict[n])
    x = np.arange(1,len(mean_mean)+1)
    
    
    plt.figure(figsize=(10,3))
    ax=plt.subplot(1,2,1)
    for i in range(len(x)):
        plt.errorbar(x[i], mean_mean[i], yerr=None, color=colors_te[i], fmt='.', 
                     capsize=3)
    
    small = dict(markersize=3)
    box1 = plt.boxplot(All_mean, flierprops=small)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors_te):
                patch.set(color=color, lw=1.7)
                patch2.set(color='k', lw=1.7)
    
    plt.xticks(labels=[], ticks=[])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    for y1, y2 in zip(All_mean[:,1], All_mean[:,2]):
        plt.plot([2,3], [y1, y2], 'darkslategray', lw=1)
    plt.annotate('Still', xy=(0.17, -0.1), xytext=(0.17, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom', 
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5),
                 arrowprops=dict(arrowstyle='-[, widthB=1.4, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    plt.annotate('Nod', xy=(0.67, -0.1), xytext=(0.67, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom',
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5), 
                 arrowprops=dict(arrowstyle='-[, widthB=3.0, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    Show_Stars(p_mean, ind, x, np.amax(All_mean, axis=0), arange_dh='diff', 
               col='black')
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
    plt.ylabel('Mean of histogram [$\\frac{\mu m^2}{s}$]')
    ax.text(-0.22, 0.9, string.ascii_lowercase[0], transform=ax.transAxes,
            size=21, weight='bold')
    
    ax=plt.subplot(1,2,2)
    for i in range(1):
        plt.errorbar(x[i], mean_std[i], yerr=None, color=colors_te[i], fmt='.', 
                     capsize=3)
    for i in range(1,len(x)):
        plt.errorbar(x[i], mean_std[i], yerr=None, label=names[i], 
                     color=colors_te[i], fmt='.', capsize=3)
    
    small = dict(markersize=3)
    box1 = plt.boxplot(All_std, flierprops=small)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors_te):
                patch.set(color=color, lw=1.7)
                patch2.set(color='k', lw=1.7)
    
    plt.xticks(labels=[], ticks=[])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    for y1, y2 in zip(All_std[:,1], All_std[:,2]):
        plt.plot([2,3], [y1, y2], 'darkslategray', lw=1)
    plt.annotate('Still', xy=(0.17, -0.1), xytext=(0.17, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom', 
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5), 
                 arrowprops=dict(arrowstyle='-[, widthB=1.4, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    plt.annotate('Nod', xy=(0.67, -0.1), xytext=(0.67, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom',
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5), 
                 arrowprops=dict(arrowstyle='-[, widthB=3.0, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    Show_Stars(p_std, ind, x, np.amax(All_std, axis=0), arange_dh='diff', 
               col='black')
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
    plt.ylabel('Std of histogram [$\\frac{\mu m^2}{s}$]')
    legend = plt.legend(loc='upper center', ncol = 2, 
                        bbox_to_anchor=(-0.2, -0.25), frameon=True)
    ax.text(-0.22, 0.9, string.ascii_lowercase[1], transform=ax.transAxes,
            size=21, weight='bold')
    
    legend.get_frame().set_linewidth(2)           
    plt.subplots_adjust(hspace=0.2, wspace=0.3)
    plt.savefig(out_dir+'ADC'+save+'.tiff', format='tiff', bbox_inches='tight', 
                dpi=200)
    plt.savefig(out_dir+'ADC'+save+'.png', format='png', bbox_inches='tight', 
                dpi=200)
    plt.show()
    
    
    # calculate mean value of ground truth scans:
    means = []
    for sub in subdir:
        if sub in ['HC_0'+str(i)+'/' for i in range(1,10)]:
            continue
        if sub in ['HC_'+str(i)+'/' for i in range(10,13)]:
            continue
        
        file = glob.glob('/data1/hannah/Registrations/'+sub+'/ADC/*OFF_STILL*.nii')[0]
        img = nib.load(file).get_fdata().astype(np.uint16)
        bm_file = glob.glob('/data1/hannah/Registrations/'+sub+'/bm*ADC*.nii')[0]
        bm = nib.load(bm_file).get_fdata().astype(np.uint16)
        
        bm_fl = bm.flatten()
        dat = img.flatten()
        img_fl = dat[bm_fl>0]
        
        means.append(np.mean(img_fl))
    
    print('Mean Values of ground truth scans:', means)
    print('Overall mean value:', np.mean(means))

    
    
    
    '''
    plt.figure(figsize=(13,3))
    ax=plt.subplot(1,3,1)
    for i in range(len(x)):
        plt.errorbar(x[i], mean_sum[i], yerr=None, color=colors_te[i], fmt='.', 
                     capsize=3)
    
    small = dict(markersize=3)
    box1 = plt.boxplot(All_sum, flierprops=small)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors_te):
                patch.set(color=color, lw=1.7)
                patch2.set(color='k', lw=1.7)
    
    plt.xticks(labels=[], ticks=[])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    for y1, y2 in zip(All_sum[:,1], All_sum[:,2]):
        plt.plot([2,3], [y1, y2], 'darkslategray', lw=1)
    plt.annotate('Still', xy=(0.17, -0.1), xytext=(0.17, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom',
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5),
                 arrowprops=dict(arrowstyle='-[, widthB=1.4, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    plt.annotate('Nod', xy=(0.67, -0.1), xytext=(0.67, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom',
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5),
                 arrowprops=dict(arrowstyle='-[, widthB=3.0, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    Show_Stars(p_sum, ind, x, np.amax(All_sum, axis=0), arange_dh='diff', 
               col='black')
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
    plt.ylabel('Mean absolute difference [$\\frac{\mu m^2}{s}$]')
    ax.text(-0.3, 0.9, string.ascii_lowercase[0], transform=ax.transAxes,
            size=21, weight='bold')
    
    ax=plt.subplot(1,3,2)
    for i in range(len(x)):
        plt.errorbar(x[i], mean_mean[i], yerr=None, color=colors_te[i], fmt='.', 
                     capsize=3)
    
    small = dict(markersize=3)
    box1 = plt.boxplot(All_mean, flierprops=small)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors_te):
                patch.set(color=color, lw=1.7)
                patch2.set(color='k', lw=1.7)
    
    plt.xticks(labels=[], ticks=[])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    for y1, y2 in zip(All_mean[:,1], All_mean[:,2]):
        plt.plot([2,3], [y1, y2], 'darkslategray', lw=1)
    plt.annotate('Still', xy=(0.17, -0.1), xytext=(0.17, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom', 
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5),
                 arrowprops=dict(arrowstyle='-[, widthB=1.4, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    plt.annotate('Nod', xy=(0.67, -0.1), xytext=(0.67, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom',
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5), 
                 arrowprops=dict(arrowstyle='-[, widthB=3.0, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    Show_Stars(p_mean, ind, x, np.amax(All_mean, axis=0), arange_dh='diff', 
               col='black')
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
    plt.ylabel('Mean of histogram [$\\frac{\mu m^2}{s}$]')
    ax.text(-0.22, 0.9, string.ascii_lowercase[1], transform=ax.transAxes,
            size=21, weight='bold')
    
    ax=plt.subplot(1,3,3)
    for i in range(1):
        plt.errorbar(x[i], mean_std[i], yerr=None, color=colors_te[i], fmt='.', 
                     capsize=3)
    for i in range(1,len(x)):
        plt.errorbar(x[i], mean_std[i], yerr=None, label=names[i], 
                     color=colors_te[i], fmt='.', capsize=3)
    
    small = dict(markersize=3)
    box1 = plt.boxplot(All_std, flierprops=small)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors_te):
                patch.set(color=color, lw=1.7)
                patch2.set(color='k', lw=1.7)
    
    plt.xticks(labels=[], ticks=[])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    for y1, y2 in zip(All_std[:,1], All_std[:,2]):
        plt.plot([2,3], [y1, y2], 'darkslategray', lw=1)
    plt.annotate('Still', xy=(0.17, -0.1), xytext=(0.17, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom', 
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5), 
                 arrowprops=dict(arrowstyle='-[, widthB=1.4, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    plt.annotate('Nod', xy=(0.67, -0.1), xytext=(0.67, -0.18), 
                 xycoords='axes fraction', fontsize=12, ha='center', va='bottom',
                 bbox=dict(boxstyle='square', fc='white', ec='grey', lw=1.5), 
                 arrowprops=dict(arrowstyle='-[, widthB=3.0, lengthB=0.7', 
                                 lw=1.5, color='grey'))
    Show_Stars(p_std, ind, x, np.amax(All_std, axis=0), arange_dh='diff', 
               col='black')
    lim = plt.gca().get_ylim()
    plt.ylim(lim[0],(lim[1]-lim[0])*1.1+lim[0])
    plt.ylabel('Std of histogram [$\\frac{\mu m^2}{s}$]')
    legend = plt.legend(loc='upper center', ncol = 2, 
                        bbox_to_anchor=(-0.8, -0.25), frameon=True)
    ax.text(-0.22, 0.9, string.ascii_lowercase[2], transform=ax.transAxes,
            size=21, weight='bold')
    
    legend.get_frame().set_linewidth(2)           
    plt.subplots_adjust(hspace=0.2, wspace=0.3)
    plt.savefig(out_dir+'ADC'+save+'.tiff', format='tiff', bbox_inches='tight', 
                dpi=200)
    plt.savefig(out_dir+'ADC'+save+'.png', format='png', bbox_inches='tight', 
                dpi=200)
    plt.show()'''
    
    
    
    
    
