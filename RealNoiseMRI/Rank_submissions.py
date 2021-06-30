'''
Script for ranking the test submissions basaed on the median rank profile for
SSIM values
'''

import numpy as np
from scipy.stats import rankdata

# define directories etc:
out_dir = 'Processed_Submissions'   # path to the folder where all intermediate
                                    # files for processing should be saved
participants = ['Team1', 'Team2']   # list of team names

test_subj = []      # list of IDs for test subjects
for i in range(1,10):
    test_subj.append('Test_0'+str(i))
for i in range(10,21):
    test_subj.append('Test_'+str(i))


def RankSSIM(task, test_subj, participants, directory):
    '''
    Calculates median rank profile for the submissions, based on SSIM values.
    Ties are treated with the method 'min': the minimum of the ranks that
    would have been assigned to all the tied values is assigned to each value.

    Parameters
    ----------
    task : str
        't1' or 't2' depending on which task to calculate the median rank
        profile for.
    test_subj : list
        List of test subjects.
    participants : list
        List of team names.
    directory : str
        Path to the directory where the results should be saved.

    Returns
    -------
    None.

    '''

    # 1)    for each test subject, rank the algorithms (optimal SSIM=1), with
    #       highest rank indicating best performance
    ranks = np.zeros((len(participants), len(test_subj)))
    for i, sub in enumerate(test_subj):
        ssim = []
        for part in participants:
            ssim.append(np.loadtxt(directory+part+'/Values_'+task+'.txt',
                                   skiprows=i+2, usecols=1, max_rows=1))
        ssim = np.array(ssim)
        ranks[:,i] = rankdata(ssim, method='min')

    # 2)    take the median rank of each alorithm across al test subjects to
    #       produce the median rank profile
    median_rank = np.median(ranks, axis=1)

    # sort the rank profile:
    ind = np.argsort(median_rank)
    median_rank_s = median_rank[ind]
    participants_s = participants[ind]

    # save the result:
    save_arr = np.array([participants_s[::-1], median_rank_s[::-1]]).T
    np.savetxt(directory+'/Ranks_'+task+'.txt', save_arr, fmt='%s',
               header='Median rank for teams specified in first column')


RankSSIM('t1', test_subj, participants, out_dir)
RankSSIM('t2', test_subj, participants, out_dir)
