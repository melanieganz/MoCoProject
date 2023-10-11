import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


def PerformWilcoxonMotion(type_metric, motion, values_metrics, out_dir, save, option=False):
    '''
    Perfrom wilcoxon signed rank test on the motion metrics.
    
    Parameters
    ----------
    type_metric : array of str
        description of metric for saving p-values.
    motion : str
        description of type of motion
    values_metrics : list of arrays
        values of the three metrics for which the test should be performed.
    option : str, optional
        If option is 'qs', then tests are performed on the quality scores. 
        The default is False.

    Returns
    -------
    p_values : array
        result of the tests.
    rej : array
        output of FDR correction - whether the hypothesis should be rejected or not.
    ind : array
        indices indicating which columns of values are compared.
    alt : array or string
        describing the alternative used for the test.

    '''
    
    text =  ['Tests RMS, median, max: ']
    p_values, stat = [], []
    
    for values, metr in zip(values_metrics, type_metric):
        #if motion in ['STILL', 'NOD']:
        if motion in ['run-01', 'run-02']:
            #ind = np.array([[0,1], [2,3], [4,5], [6,7], [8,9], [10,11]]) #det sidste svarer til diff, så det skal tilbage når diff er klar
            ind = np.array([[0,1], [2,3], [4,5], [6,7], [8,9]])
            alt = 'two-sided'
                
        #elif motion == 'SHAKE':
        elif motion == 'run-03':
            ind = np.array([[0,1]])
            alt = 'two-sided'
        
        for index in ind:
            print('motion', motion)
            print('index', index)
            print('values', values)
            tmp, p = wilcoxon(values[index[0]], values[index[1]], alternative=alt)
            text.append(metr+' '+str(index[0])+' '+str(index[1])+' '+alt+': p = '+str(p))
            p_values.append(p)
            stat.append(tmp)
    
    text = np.array(text)
    if not os.path.exists(out_dir+'MotionMetrics_'+motion+'/'):
        os.makedirs(out_dir+'MotionMetrics_'+motion+'/')
    np.savetxt(out_dir+'MotionMetrics_'+motion+'/ParametricTests'+save+'.txt', text, fmt='%s')
    
    # effect size:
    text =  ['Tests RMS, median, max: ']
    stat = np.array(stat)
    n = []
    # sequ = ['MPR', 'TSE', 'TIRM', 'FLAIR', 'T2STAR', 'DWI']
    sequ = ['mprage', 't2tse', 't1tirm', 'flair', 't2star']
    for j in [1,2,3]:
        #if motion in ['STILL', 'NOD']:
        if motion in ['run-01', 'run-02']:
            for i in range(0,len(sequ)):
                n.append(len(values_metrics[0][2*i]))
        #elif motion == 'SHAKE':
        elif motion == 'run-03':
            n.append(len(values_metrics[0][0]))
    n = np.array(n)
    mu = n*(n+1)/4
    s = np.sqrt(n*(n+1)*(2*n+1)/24)
    z = (stat-mu)/s
    ES = np.abs(z)/np.sqrt(n)

    for j in range(0,len(n)):
        text.append(str(j)+' '+sequ[j%len(sequ)]+' n:'+str(n[j])+' ES: '+str(ES[j]))
    
    print(n,ES)
    np.savetxt(out_dir+'MotionMetrics_'+motion+'/EffectSize'+save+'.txt', text, fmt='%s')
       
    return p_values, ind, alt, ES



def PerformWilcoxonAllImg(type_metric, values, sequ, out_dir, save, option=False):
    '''
    Performs wilcoxon signed rank test for metrics claculated on all 
    propsectively corrected and uncorrected images.

    Parameters
    ----------
    type_metric : str
        description of metric for saving p-values.
    values : array
        values of the metric for which the test should be performed.
    sequ : str
        description / abbr. of the sequence.
    option : str, optional
        If option is 'qs', then tests are performed on the quality scores.
        If option is 'diff', then indices are adjusted accordingly.
        The default is False.

    Returns
    -------
    p_values : array
        result of the tests.
    rej : array
        output of FDR correction - whether the hypothesis should be rejected or not.
    ind : array
        indices indicating which columns of values are compared.
    alt : array or string
        describing the alternative used for the test.

    '''

    text =  ['Tests '+type_metric+': ']
    p_values = []

    #if sequ == 'T1_MPR':
    if sequ == 'mprage':
        ind = np.array([[0,1], [2,3], [4,5], [2,5],[2,4], [3,5], [6,7], [8,9], [6,9], [6,8], [7,9]])
    

    #elif sequ in ['T2STAR', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
    elif sequ in ['t2star']:
        ind = np.array([[0,1],[2,3]])

        #if option == 'diff':
        #    ind = np.array([[0,1], [1,2], [0,2]])
    else:
        ind = np.array([[0,1],[2,3], [4,5], [2,5], [2,4], [3,5]])
        
    altern = 'two-sided' 
    
    #if type_metric=='QS' and sequ == 'T2STAR':
    if type_metric=='QS' and sequ == 't2star':
        p_values = [0,0]  
    else:
        for index in ind:        
            # sort out nans:
            values_1 = values[:,index[0]]
            values_2 = values[:,index[1]]
            values_1 = values_1[np.logical_not(np.isnan(values_1))]
            values_2 = values_2[np.logical_not(np.isnan(values_2))]
            
            if len(values_1) != len(values_2):
                print('ERROR: arrays for '+type_metric+' contain different number of NaN elements')
            
              
            if option == 'qs':
                tmp, p = wilcoxon(values_1, values_2, alternative=altern, zero_method='pratt') # not used anymore!!!
            else:
                tmp, p = wilcoxon(values_1, values_2, alternative=altern)
            text.append(str(index[0])+' '+str(index[1])+' '+altern+': p = '+str(p))
            p_values.append(p)

    #correct for multiple comparisons:
    rej, p_values_cor, _, __ = multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    text.append('Corrected p-values')
    text.append(p_values_cor.astype(str))

    text = np.array(text, dtype=object)
    #if option == 'diff':
    #    tmp = out_dir+'ParametricTests/'+sequ+'/Tests_diff_'+type_metric+save+'.txt'
    #    if not os.path.exists(out_dir+'ParametricTests/'+sequ+'/'):
    #        os.makedirs(out_dir+'ParametricTests/'+sequ+'/')
    #    np.savetxt(tmp, text, fmt='%s')
    #else:
    #    if not os.path.exists(out_dir+'ParametricTests/'+sequ+'/'):
    #        os.makedirs(out_dir+'ParametricTests/'+sequ+'/')
    #    np.savetxt(out_dir+'ParametricTests/'+sequ+'/Tests_'+type_metric+save+'.txt', text, fmt='%s')

    return p_values_cor, rej, ind, altern



def PerformWilcoxonRetro(type_metric, values, sequ, out_dir, save, option=False):
    '''
    Performs wilcoxon signed rank test for metrics claculated on all 
    MPRAGE and FLAIR scans from retrospective reconstruction pipeline.

    Parameters
    ----------
    type_metric : str
        description of metric for saving p-values.
    values : array
        values of the metric for which the test should be performed.
    sequ : str
        description / abbr. of the sequence.
    option : str, optional
        If option is 'qs', then tests are performed on the quality scores.
        The default is False.

    Returns
    -------
    p_values : array
        result of the tests.
    rej : array
        output of FDR correction - whether the hypothesis should be rejected or not.
    ind : array
        indices indicating which columns of values are compared.
    alt : array or string
        describing the alternative used for the test.

    '''

    text =  ['Tests '+type_metric+': ']
    p_values = []

    #if sequ == 'T1_MPR':
    if sequ == 'mprage':
        ind = np.array([[0,1], [1,2], [0,2], [3,4], [4,5], [3,5], [6,7], [7,8], [6,8], [3,8], [3,6], [4,7], [5,8], [9,10], [10,11], [9,11], [12,13], [13,14], [12,14], [9,14], [9,12], [10,13], [11,14]])


    #elif sequ in ['T2STAR', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
    elif sequ in ['t2star']:
        ind = np.array([[0,1],[2,3]])

    else:
        ind = np.array([[0,1],[1,2],[0,2],[3,4],[4,5],[3,5],[6,7],[7,8],[6,8],[3,8],[3,6], [4,7], [5,8]])

    for index in ind:
        altern = 'two-sided'
        if option == 'qs':
            tmp, p = wilcoxon(values[:,index[0]], values[:,index[1]], alternative=altern, zero_method='pratt') #not used anymore
        else:
            tmp, p = wilcoxon(values[:,index[0]], values[:,index[1]], alternative=altern)
        text.append(str(index[0])+' '+str(index[1])+' '+altern+': p = '+str(p))
        p_values.append(p)


    #correct for multiple comparisons:
    rej, p_values_cor, _, __ = multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    text.append('Corrected p-values')
    text.append(p_values_cor)

    text = np.array(text)
    np.savetxt(out_dir+'ParametricTests/'+sequ+'/Tests_'+type_metric+save+'.txt', text, fmt='%s')

    return p_values_cor, rej, ind, altern
