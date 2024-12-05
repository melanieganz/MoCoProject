import os
import pandas as pd
import numpy as np
import string
from statsmodels.stats.multitest import multipletests
from utils import Show_Stars, SortFiles, MakeBoxplot, DrawLines2, DrawLines
from statistical_tests import PerformWilcoxonMotion, PerformWilcoxonAllImg, PerformWilcoxonAllImg_custom_indices
from functools import reduce
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
import matplotlib


sns.set_style("whitegrid")
#plt.style.use('seaborn-whitegrid')
matplotlib.rc('axes', edgecolor='black')



root_dir = os.environ.get("MOCO_DATASET_PATH")
out_dir = os.path.join(root_dir, "derivatives/results/plots")
in_dir_motion_estimates = os.path.join(root_dir, "derivatives/results/Motion_Estimates/")
in_dir_metric_results = os.path.join(root_dir, "derivatives/results/metricsresults/")
in_dir_observer_scores =  os.path.join(root_dir, "derivatives/observer_scores")
out_dir_metrics = os.path.join(root_dir, "derivatives/results/Metrics_Results/Comparison/")

save = '_2022_05_27'

SMALL_MARKER = dict(markersize=3)



sequs = ['mprage', 'flair', 't2tse', 't1tirm', 't2star']

bids_seq_to_observer_score_map = {
    "mprage": "acq-mprage_T1w.tsv",
    "flair":  "acq-flair_FLAIR.tsv",
    "t2tse":  "acq-t2tse_T2w.tsv",
    "t1tirm": "acq-t1tirm_T1w.tsv",
}

dwi_sequence = "dwi"
dwi_filename = "dwi.tsv"

GT_STILL_IMAGE_SUBSEQUENCE = "pmcoff_rec-wore_run-01"


def read_scores():
    seq_to_df = {}
    seq_to_df_raw = {}
    for seq, file_name in bids_seq_to_observer_score_map.items():
        file_path = os.path.join(in_dir_observer_scores, file_name)
        scores_df : pd.DataFrame = pd.read_csv(file_path, delimiter='\t')

        # combine scores into one score
        score_name = "score"
        scores_df[score_name] = scores_df.apply(lambda row: (row["radiographer_1_score"] + row["radiographer_2_score"]
                                            + 2 * row["neuroradiologist_score"]) / 4, axis=1)
        scores_df.drop(["radiographer_1_score", "radiographer_2_score", "neuroradiologist_score"], axis=1, inplace=True)
        
        scores_df["sub_sequence"] = scores_df["sequence"].apply(func= lambda seq: seq[seq.find("pmc"):seq.find("run-") + 6])
        scores_df.drop("sequence", inplace=True, axis=1)
        scores_df.rename({"participant_id": "subject"}, inplace=True, axis=1)

        seq_to_df_raw[seq] = scores_df.copy(deep=True)

        scores_df.rename({score_name : "score_" + seq}, inplace=True, axis=1)
        score_name = "score_" + seq

        scores_df = scores_df[["sub_sequence", score_name]].groupby("sub_sequence").agg("mean")
        
        seq_to_df[seq] = scores_df

    def sort_ignoring_pmc(index):
        return [i[i.find("rec"):] for i in index]


    joined_score_df : pd.DataFrame= reduce(lambda df1, df2: df1.merge(df2, "inner", on="sub_sequence"), seq_to_df.values())
    joined_score_df.sort_index(key=sort_ignoring_pmc, inplace=True)

    return joined_score_df, seq_to_df_raw

def read_dwi_rank():
    

    file_path = os.path.join(in_dir_observer_scores, dwi_filename)
    rank_df : pd.DataFrame = pd.read_csv(file_path, delimiter='\t')

    rank_df["rank"] = rank_df["neuroradiologist_rank"]
    rank_df.drop("neuroradiologist_rank", axis=1, inplace=True)
    
    rank_df["sub_sequence"] = rank_df["sequence"].apply(func= lambda seq: seq[seq.find("pmc"):seq.find("run-") + 6])
    rank_df.drop("sequence", inplace=True, axis=1)

    rank_df.rename({"participant_id": "subject"}, inplace=True, axis=1)

    return rank_df

def read_image_metrics():
    seq_to_img_metrics = {}
    for seq in sequs:
        img_metrics_df : pd.DataFrame = pd.DataFrame()
        for sub in [f"sub-{i:02d}" for i in range(1, 23)]:
            metrics = {}

            sub_folder_path = os.path.join(in_dir_metric_results, sub)
            all_files = os.listdir(sub_folder_path)
            seq_files = [file for file in all_files if seq in file]
            if len(seq_files) == 0:
                continue
            if len(seq_files) > 1:
                raise Exception("Expected exactly one metrics file for each sequence")
            
            seq_file = os.path.join(sub_folder_path,seq_files[0])
        
            sub_sequences = np.loadtxt(seq_file, unpack=True, dtype=str, usecols=0, skiprows=1)
            ssim, psnr, tg = np.loadtxt(seq_file, unpack=True, usecols=(1,2,3), skiprows=1)

            metrics["subject"] = sub
            metrics["sub_sequence"] = sub_sequences
            metrics["ssim"] = ssim
            metrics["psnr"] = psnr
            metrics["tg"] = tg


            img_metrics_df = pd.concat([img_metrics_df, pd.DataFrame.from_dict(metrics)])

        # normalise tg values by dividing with mean value of off still
        mean_still_tg = img_metrics_df.groupby("sub_sequence").agg({"tg" : "mean"}).loc[GT_STILL_IMAGE_SUBSEQUENCE]
        img_metrics_df["tg"] = img_metrics_df["tg"].apply(lambda tg: tg / mean_still_tg)

        
        seq_to_img_metrics[seq] = img_metrics_df
    return seq_to_img_metrics

# contains df of all subjects and subsequences with metric scores, contained in a map by sequence (MPR, FLAIR, etc.)
def calculate_p_values_observer_scores(relevant_subsequences, seq_to_observer_scores_df, indices=[[0, 1]]):
    seq_to_observer_scores_p_values = {}

    for seq, observer_scores_df in seq_to_observer_scores_df.items():
        grouped_df = observer_scores_df.groupby("sub_sequence").agg(list)

        relevant_subsequence_scores = grouped_df.loc[relevant_subsequences]
        scores_arr = np.array(relevant_subsequence_scores["score"].to_list()).T
        p_qs, rej_qs, ind_p, alt = PerformWilcoxonAllImg_custom_indices('QS', scores_arr, seq, out_dir_metrics, save, ind=indices)

        seq_to_observer_scores_p_values[seq] = {tuple(np.array(relevant_subsequences)[indices[i]]): p_q for i, p_q in enumerate(p_qs)}
    return seq_to_observer_scores_p_values


def calculate_p_values_observer_rank(relevant_subsequences, dwi_rank_df):

    relevant_subsequence_ranks = dwi_rank_df.groupby("sub_sequence").agg(list).loc[relevant_subsequences]
    ranks_arr = np.array(relevant_subsequence_ranks["rank"].to_list()).T

    p_qs, rej_qs, ind_p, alt = PerformWilcoxonAllImg_custom_indices('QS', ranks_arr, "dwi", out_dir_metrics, save, ind=[[0, 1]])

    return dict(zip([tuple(relevant_subsequences)] * len(p_qs), p_qs))

def calculate_p_values_image_metrics(relevant_subsequences, seq_to_img_metrics, indices=[[0, 1]] ):
    metric_to_seq_to_p_values = {}

    for metric in ( "ssim", "psnr", "tg"):
        seq_to_p_values = {}
        for seq, image_metric_df in seq_to_img_metrics.items():
            grouped_df = image_metric_df.groupby("sub_sequence").agg(list)

            relevant_subsequence_metrics = grouped_df.loc[relevant_subsequences]
            metric_arr = np.array(relevant_subsequence_metrics[metric].to_list()).T
            if len(indices) == 0:
                p_qs = []
            else:
                p_qs, rej_qs, ind_p, alt = PerformWilcoxonAllImg_custom_indices(metric.upper(), metric_arr, seq, out_dir_metrics, save, ind=indices)
            
            seq_to_p_values[seq] = {tuple(np.array(relevant_subsequences)[indices[i]]): p_q for i, p_q in enumerate(p_qs)}

        metric_to_seq_to_p_values[metric] = seq_to_p_values
    return metric_to_seq_to_p_values


def draw_boxplot_1(seq_to_img_metrics, seq_to_df, dwi_rank_df, relevant_subsequences):
    seq_to_observer_scores_p_values = calculate_p_values_observer_scores(relevant_subsequences, seq_to_df)

    image_metrics_p_values = calculate_p_values_image_metrics(relevant_subsequences, seq_to_img_metrics)

    # removing reacquisition from subsequences since it was not part of DWI
    dwi_subsequences = ["_".join([subs.split("_")[0]] + [subs.split("_")[2]]) for subs in relevant_subsequences]
    dwi_p_values = calculate_p_values_observer_rank(dwi_subsequences, dwi_rank_df)

    labels = ['Tenengrad', 'Observer Scores', 'Quality Rank']
    colors = ['tab:orange', 'tab:blue', 'tab:orange', 'tab:blue', 'tab:orange', 
              'tab:blue', 'tab:orange', 'tab:blue', 'tab:orange', 'tab:blue', 
              'tab:orange', 'tab:blue']
    fig_labels = ['without PMC', 'with PMC']
    x_axis_range = np.arange(1,13)
    x_axis_even_range = np.arange(0, 13, 2)
    x_axis_odd_range = np.arange(1, 12, 2)


    def draw_tenengrad_image_metric():
        relevant_metrics_for_sequence = []
        means = []
        for seq in sequs:
            df : pd.DataFrame = seq_to_img_metrics[seq]
            df = df[df["sub_sequence"].isin(relevant_subsequences)]
            relevant_metrics_for_sequence += [df[df["sub_sequence"] == subseq]["tg"].to_numpy() for subseq in relevant_subsequences]
            means += df.groupby("sub_sequence").agg({"tg": "mean"})["tg"].to_list()


        ax = plt.subplot2grid((2,6), (0,0), colspan=6)
        box1 = plt.boxplot(relevant_metrics_for_sequence, flierprops=SMALL_MARKER)
        for j in range(len(means)):
            plt.errorbar(x_axis_range[j], means[j], yerr=None, color=colors[j], fmt='.', capsize=3)
        ticklabels = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_STIR', 'T2*']
        ticks = [1.5, 3.5, 5.5, 7.5, 9.5] # when DIFF added, add 11.5 to list
        plt.xticks(labels=ticklabels, ticks=ticks, fontsize=14)

        for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
            patch.set(color=color, lw=1.7)
            patch2.set(color='k', lw=1.7)
        
        def draw_p_values(relevant_subsequences, relevant_metric_for_sequence):

            def draw_stars(m_still, p_values):
                maxi = []
                for v in m_still:
                    maxi.append(np.amax(v))

                indices = [[0,1], [2,3], [4,5], [6,7], [8,9]]
                Show_Stars(np.array(p_values), indices[0:len(p_values)], x_axis_range, maxi,
                            col='black')

            def draw_lines(m_still, p_values):
                DrawLines2(x_axis_even_range[0:len(p_values)],x_axis_odd_range[0:len(p_values)], x_axis_odd_range[0:len(p_values)],
                        x_axis_even_range[1:len(p_values) + 1], m_still, lw=0.7, col='darkslategray')
            
            p_values = [x for p_vals in image_metrics_p_values["tg"].values() for x in list(p_vals.values())]
            draw_stars(relevant_metric_for_sequence, p_values)
            draw_lines(relevant_metric_for_sequence, p_values)

        draw_p_values(relevant_subsequences, relevant_metrics_for_sequence)
        
        plt.ylabel("Tenengrad", fontsize=15)
        plt.yticks(fontsize=13)
        plt.tick_params('both', length=0)
        plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    


    
    def draw_observer_scores(): 

        relevant_scores_for_sequence = []
        means = []
        for seq in bids_seq_to_observer_score_map.keys():
            df : pd.DataFrame = seq_to_df[seq]
            df = df[df["sub_sequence"].isin(relevant_subsequences)]
            relevant_scores_for_sequence += [df[df["sub_sequence"] == subseq]["score"].to_numpy() for subseq in relevant_subsequences]
            means += df.groupby("sub_sequence").agg({"score": "mean"})["score"].to_list()

        ax = plt.subplot2grid((2,6), (1,0), colspan=4)
        box1 = plt.boxplot(relevant_scores_for_sequence, flierprops=SMALL_MARKER)
        for j in range(len(means)):
            plt.errorbar(x_axis_range[j], means[j], yerr=None, color=colors[j], fmt='.', capsize=3)
        for j in range(0,2):
            plt.errorbar(x_axis_range[j], means[j], yerr=None, color=colors[j], fmt='.', capsize=3, label=fig_labels[j])
        ticklabels = ['T1_MPR', 'T2_FLAIR', 'T2_TSE', 'T1_STIR']
        ticks = [1.5, 3.5, 5.5, 7.5]
        plt.xticks(labels=ticklabels, ticks=ticks, fontsize=14)

        for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
            patch.set(color=color, lw=1.7)
            patch2.set(color='k', lw=1.7)


        def draw_p_values(relevant_sequences, relevant_scores_for_sequence):

            def draw_stars(m_still, p_values):
                maxi = []
                for v in m_still:
                    maxi.append(np.amax(v))

                indices = [[0,1], [2,3], [4,5], [6,7]]
                Show_Stars(np.array(p_values), indices[0:len(p_values)], x_axis_range[0:8], maxi[0:8],
                            col='black')

            def draw_lines(m_still, p_values):
                # a, b, c, d
                DrawLines2(x_axis_even_range[0:len(p_values)],x_axis_odd_range[0:len(p_values)], x_axis_odd_range[0:len(p_values)],
                        x_axis_even_range[1:len(p_values) + 1], m_still, lw=0.7, col='darkslategray')
            
            p_values = [p for seq in relevant_sequences for p in seq_to_observer_scores_p_values[seq].values() ]
            draw_stars(relevant_scores_for_sequence, p_values)
            draw_lines(relevant_scores_for_sequence, p_values)
        
        

        draw_p_values(bids_seq_to_observer_score_map.keys(), relevant_scores_for_sequence)
        
        plt.ylabel("Observer Scores", fontsize=15)
        plt.yticks(fontsize=13)
        plt.tick_params('both', length=0)
        plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))


        plt.yticks(ticks=[2.5, 3, 3.5, 4, 4.5, 5])
        ax.text(-0.2, 0.95, string.ascii_lowercase[1],
                    transform=ax.transAxes, size=24, weight='bold')
    
    def draw_rank(relevant_subsequences):

        dwi_rank_df_relevant = dwi_rank_df[dwi_rank_df["sub_sequence"].isin(relevant_subsequences)]
        relevant_ranks_for_sequence = [dwi_rank_df_relevant[dwi_rank_df_relevant["sub_sequence"] == subseq]["rank"].to_numpy() for subseq in relevant_subsequences]
        means = dwi_rank_df_relevant.groupby("sub_sequence").agg({"rank": "mean"})["rank"].to_list()

        ax = plt.subplot2grid((2,6), (1,5), colspan=1)
        box1 = plt.boxplot(relevant_ranks_for_sequence, flierprops=SMALL_MARKER)
        for j in range(len(means)):
            plt.errorbar(x_axis_range[j], means[j], yerr=None, color=colors[j], fmt='.', capsize=3)
        ticklabels = ['DWI']
        ticks = [1.5]
        plt.xticks(labels=ticklabels, ticks=ticks, fontsize=14)
        plt.yticks([2,3,4])

        def draw_p_values(m_still, p_still):
            DrawLines2(x_axis_even_range[0:len(p_still)],x_axis_odd_range[0:len(p_still)], x_axis_odd_range[0:len(p_still)],
                        x_axis_even_range[1:len(p_still) + 1], m_still, lw=0.7, col='darkslategray')
        
        draw_p_values(relevant_ranks_for_sequence, list(dwi_p_values.values()))

        for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
            patch.set(color=color, lw=1.7)
            patch2.set(color='k', lw=1.7)
        
        
        plt.ylabel("Quality Rank", fontsize=15)
        plt.yticks(fontsize=13)
        plt.tick_params('both', length=0)
        plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.text(-1, 0.95, string.ascii_lowercase[2],
                    transform=ax.transAxes, size=24, weight='bold')

    plt.figure(figsize=(10,9))

    draw_tenengrad_image_metric()
    draw_observer_scores()
    draw_rank(dwi_subsequences)

    legend = plt.legend( loc='lower left', ncol=2,
                                bbox_to_anchor=(0.4, -0.3), fontsize=14,
                                frameon=True)
    plt.tight_layout()
    legend.get_frame().set_linewidth(2)
    plt.subplots_adjust(hspace=0.2, wspace=0.3)
    plt.show()


# still plot
relevant_subsequences_plot_1 = ["pmcoff_rec-wore_run-01", "pmcon_rec-wore_run-01"]
seq_to_img_metrics = read_image_metrics()
_, seq_to_observer_scores_df = read_scores()
dwi_rank_df = read_dwi_rank()

draw_boxplot_1(seq_to_img_metrics, seq_to_observer_scores_df, dwi_rank_df, relevant_subsequences_plot_1)

# nodding plot

relevant_subsequences_plot_2 = {"mprage": ["pmcoff_rec-wore_run-03", "pmcon_rec-wore_run-03", "pmcoff_rec-wre_run-03", "pmcon_rec-wre_run-03",
                                           "pmcoff_rec-wore_run-02", "pmcon_rec-wore_run-02", "pmcoff_rec-wre_run-02", "pmcon_rec-wre_run-02"],
                                          
                                "flair" :["pmcoff_rec-wore_run-02", "pmcon_rec-wore_run-02", "pmcoff_rec-wre_run-02", "pmcon_rec-wre_run-02"]
                                }


# common box plot variables
labels = ['SSIM', 'PSNR', 'Tenengrad', 'Observer Scores', "DWI Rank"]
color_dict = {
                  'no MoCo no REAC':'tab:orange',
                  'prospective MoCo no REAC': 'tab:blue', 
                  'no MoCo with REAC':'tab:cyan',
                  'prospective MoCo with REAC': 'tab:green', 
                  }
colors = list(color_dict.values()) * 4


def draw_boxplot_2(relevant_sequences, sequence_to_relevant_subsequences):
    
    x = np.concatenate((np.arange(1,9), np.arange(10,14)), axis=None)
    a = [0,2,4,6] 
    b = [1,3,5,7]
    c = [1,3,5,7]
    d = [2,4,6,8]
    a_ = [0,2]
    b_ = [1,3]
    c_ = [10,12]
    d_ = [11,13]

    plt.figure(figsize=(11,8))

    for i, metric in enumerate(["ssim", "psnr", "tg"]):

        for seq, positions, x_axis_for_lines in zip(relevant_sequences, [None, range(10,14)], [(a,b,c,d), (a_,b_,c_,d_)]):
            relevant_subsequences = sequence_to_relevant_subsequences[seq]
            indices = [[0, 1], [2, 3], [0, 3], [0, 2], [1, 3],
                       [4, 5], [6, 7], [4, 7], [4, 6], [5, 7]
                      ]
            if seq == "flair":
                indices = indices[0:5]
             
            p_values = calculate_p_values_image_metrics(relevant_subsequences, {k:v for k,v in seq_to_img_metrics.items() if k == seq},
                                                                   indices=indices)[metric][seq]

            df : pd.DataFrame = seq_to_img_metrics[seq]
            df = df[df["sub_sequence"].isin(relevant_subsequences)]

            relevant_metrics_for_sequence = [df[df["sub_sequence"] == subseq][metric].to_numpy() for subseq in relevant_subsequences]
            means = np.mean(relevant_metrics_for_sequence, axis=1)
            
            ax = plt.subplot(2,2,i+1)
            draw_subgraph(colors, x, i, positions, 
                               x_axis_for_lines, indices, np.array(list(p_values.values())),
                               relevant_metrics_for_sequence, means, ax)
            lim = plt.gca().get_ylim()
            plt.ylim(lim[0],(lim[1]-lim[0])*1.01+lim[0])


            
        
        subplot_settings_plot_2(labels, i, ax)

    # observer scores
    for seq, positions, x_axis_for_lines in zip(relevant_sequences, [None, range(10,14)], [(a,b,c,d), (a_,b_,c_,d_)]):
        relevant_subsequences = sequence_to_relevant_subsequences[seq]
        indices = [[0, 1], [2, 3], [0, 3], [0, 2], [1, 3],
                    [4, 5], [6, 7], [4, 7], [4, 6], [5, 7]
                    ]
        if seq == "flair":
            indices = indices[0:5]
        
        p_values, relevant_metrics_for_sequence, means = \
            get_observer_score_values_for_plotting(seq, relevant_subsequences, indices)
        
        ax = plt.subplot(2,2,4)
        draw_subgraph(colors, x, 3, positions, x_axis_for_lines,
                            indices, np.array(list(p_values.values())), relevant_metrics_for_sequence, means, ax)
        lim = plt.gca().get_ylim()
        plt.ylim(lim[0],(lim[1]-lim[0])*1.01+lim[0])
    
    subplot_settings_plot_2(labels, 3, ax)
    plt.yticks(ticks=[1, 2, 3, 4, 5])
    

    legend = plt.legend( loc='lower left', ncol=2, 
                    bbox_to_anchor=(0.25, -0.4), fontsize=12, 
                    frameon=True)
    legend.get_frame().set_linewidth(2)        

    plt.subplots_adjust(hspace=0.2, wspace=0.3)
    plt.show()

def get_observer_score_values_for_plotting(seq, relevant_subsequences, indices):
    p_values = calculate_p_values_observer_scores(relevant_subsequences, 
                                                                    {k:v for k,v in seq_to_observer_scores_df.items() if k == seq},
                                                                    indices=indices)[seq]

    df : pd.DataFrame = seq_to_observer_scores_df[seq]
    df = df[df["sub_sequence"].isin(relevant_subsequences)]
        # this uses the order of relevant_subsequences, which is important so it matches the graph
    relevant_metrics_for_sequence = [df[df["sub_sequence"] == subseq]["score"].to_numpy() for subseq in relevant_subsequences]
    means = np.mean(relevant_metrics_for_sequence, axis=1)
    return p_values,relevant_metrics_for_sequence,means

def subplot_settings_plot_2(labels, i, ax):
    plt.ylabel(labels[i], fontsize=15)

    ticklabels = ['T1_MPR\nSHAKE', 'T1_MPR\nNOD', 'T2_FLAIR\nNOD']
    ticks = [2.5, 6.5, 11.5]
    plt.xticks(labels=ticklabels, ticks=ticks, fontsize=12)
    ax.text(-0.18, 0.9, string.ascii_lowercase[i], transform=ax.transAxes,
                size=21, weight='bold')
    plt.yticks(fontsize=13)
    plt.tick_params('both', length=0)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

def draw_subgraph(colors, x, i, positions, x_axis_for_lines, indices, p_values, relevant_metrics_for_sequence, means, ax):

    box1 = plt.boxplot(relevant_metrics_for_sequence, positions=positions, flierprops=SMALL_MARKER, widths=0.5)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
        patch.set(color=color, lw=1.5)
        patch2.set(color='k', lw=1.5)

    for j in range(len(means)):
        plt.errorbar((x if positions is None else positions)[j], means[j], yerr=None, color=colors[j], fmt='.', 
                            capsize=3)

    def draw_p_values():
        maxi = []
        for v in np.array(relevant_metrics_for_sequence):
            maxi.append(np.amax(v))

        Show_Stars(p_values, np.array(indices), x if positions is None else positions, maxi, 
                        arange_dh='PAPER', col='black')
            
    draw_p_values()
    DrawLines(*x_axis_for_lines, np.array(relevant_metrics_for_sequence).T, col='darkslategray')
    return ax
    
relevant_sequences_plot_2 = ['mprage', 'flair']

draw_boxplot_2(relevant_sequences_plot_2, relevant_subsequences_plot_2)



def draw_boxplot_3(relevant_sequences, seq_to_relevant_subsequences):
    x_axis_range = np.arange(0,16)
    x_axis_even_range = np.arange(0, 17, 2)
    x_axis_odd_range = np.arange(1, 16, 2)

    plt.figure(figsize=(15,8))

    indices = [[[0, 1], [2, 3], [0, 3], [0, 2], [1, 3]],
               [[0, 1], [2, 3], [0, 3], [0, 2], [1, 3]],
               []
               ]


    metric_sequences = [x for x in relevant_sequences if not x == "TRACEW"]

    for i, metric in enumerate(["ssim", "psnr", "tg"]):
        # collected for each sequence
        relevant_metrics_for_sequence = []
        lims = []
        
        ax = plt.subplot2grid((2,4), (i//2,i%2*2), colspan=2)

        for in_plot_offset, seq, ind in zip(range(1, 12, 5), metric_sequences, indices):
            relevant_subsequences = seq_to_relevant_subsequences[seq]
            
            p_values = list(calculate_p_values_image_metrics(relevant_subsequences, 
                                                {k:v for k,v in seq_to_img_metrics.items() if k == seq},
                                                indices=ind)[metric][seq].values())
            df : pd.DataFrame = seq_to_img_metrics[seq]
            df = df[df["sub_sequence"].isin(relevant_subsequences)]

            relevant_metrics_for_sequence = [df[df["sub_sequence"] == subseq][metric].to_numpy() for subseq in relevant_subsequences]
            means = np.mean(relevant_metrics_for_sequence, axis=1)

            

            positions = range(in_plot_offset, in_plot_offset + len(relevant_subsequences))
            a = x_axis_even_range[:len(relevant_subsequences) // 2]
            b = x_axis_odd_range[:len(relevant_subsequences) // 2]
            c = [positions[p] for p in a]
            d = [positions[p] for p in b]
            x_axis_for_lines = (a,b,c,d)
            ax = draw_subgraph(colors, x_axis_range, i, positions, 
                          x_axis_for_lines, 
                          ind, np.array(p_values),
                          relevant_metrics_for_sequence, means, ax)
            
            lims.append(plt.gca().get_ylim())
        
        subplot_settings_plot_3(labels, i, ax)
        lims = np.array(lims)
        lim_min = np.amin(lims[:,0])
        lim_max = np.amax(lims[:,1])
        plt.ylim(lim_min,(lim_max-lim_min)*1.01+lim_min)

    def draw_quality_score():
        ax = plt.subplot2grid((2,6), (1,3), colspan=2)
        qs_relevant_sequences = [s for s in relevant_sequences if not s == "TRACEW" and not s == "t2star"]
        for in_plot_offset, seq, ind in zip(range(1, 12, 5), qs_relevant_sequences, indices):
            relevant_subsequences = seq_to_relevant_subsequences[seq]

            p_values, relevant_metrics_for_sequence, means = \
                get_observer_score_values_for_plotting(seq, relevant_subsequences, ind)
            
            positions = range(in_plot_offset, in_plot_offset + len(relevant_subsequences))
            a = x_axis_even_range[:len(relevant_subsequences) // 2]
            b = x_axis_odd_range[:len(relevant_subsequences) // 2]
            c = [positions[p] for p in a]
            d = [positions[p] for p in b]
            x_axis_for_lines = (a,b,c,d)
            draw_subgraph(colors, x_axis_range, 3, positions, x_axis_for_lines,
                                ind, np.array(list(p_values.values())), relevant_metrics_for_sequence, means, ax)
            
            
        
        subplot_settings_plot_3(labels, 3, ax, num_seqs=len(qs_relevant_sequences))

        plt.yticks([1,2,3,4,5])
        lim = plt.gca().get_ylim()
        plt.ylim(lim[0],(lim[1]-lim[0])*1.01+lim[0])

    draw_quality_score()
            
    def draw_quality_rank():
        relevant_subsequences = seq_to_relevant_subsequences["TRACEW"]
        relevant_subsequences = ["_".join([subs.split("_")[0]] + [subs.split("_")[2]]) for subs in relevant_subsequences]
        
        ax = plt.subplot2grid((2,6), (1,5), colspan=1)  

        df: pd.DataFrame = dwi_rank_df[dwi_rank_df["sub_sequence"].isin(relevant_subsequences)]
        p_values = calculate_p_values_observer_rank(relevant_subsequences, dwi_rank_df)
        relevant_metrics_for_sequence = [df[df["sub_sequence"] == subseq]["rank"].to_numpy() for subseq in relevant_subsequences]
        means = np.mean(relevant_metrics_for_sequence, axis=1)



        in_plot_offset = 4
        positions =  np.arange(1, 1 + len(relevant_subsequences))
        a = x_axis_even_range[:len(relevant_subsequences) // 2]
        b = x_axis_odd_range[:len(relevant_subsequences) // 2]
        c = [positions[p] for p in a]
        d = [positions[p] for p in b]
        x_axis_for_lines = (a,b,c,d)
              
        draw_subgraph(colors, x_axis_range, 4, positions, x_axis_for_lines,
                            [[0,1]], np.array(list(p_values.values())), relevant_metrics_for_sequence, means, ax)

        ticklabels = ['DWI']
        ticks = [1.5]
        plt.xticks(labels=ticklabels, ticks=ticks, fontsize=12)
        plt.yticks(ticks=[1, 2, 3, 4])

    draw_quality_rank()

    
    legend = plt.legend( loc='lower left', ncol=2, 
                    bbox_to_anchor=(0.25, -0.4), fontsize=12, 
                    frameon=True)
    legend.get_frame().set_linewidth(2)        

    plt.subplots_adjust(hspace=0.2, wspace=0.6)

    plt.show()

def subplot_settings_plot_3(labels, i, ax, num_seqs=3):
    plt.ylabel(labels[i], fontsize=15)

    ticklabels = ['T2_TSE', 'T1_STIR', 'T2*'][:num_seqs]
    ticklabels = ['T2_TSE', 'T1_STIR', 'T2*'][:num_seqs]
    ticks = [2.5, 7.5, 11.5][:num_seqs]
    plt.xticks(labels=ticklabels, ticks=ticks, fontsize=12)
    ax.text(-0.18, 0.9, string.ascii_lowercase[i], transform=ax.transAxes,
                size=21, weight='bold')
    plt.yticks(fontsize=13)
    plt.tick_params('both', length=0)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

relevant_sequences_plot_3 = ['t2tse', 't1tirm', 't2star']

seq_to_relevant_subsequences = {seq:["pmcoff_rec-wore_run-02", "pmcon_rec-wore_run-02",
                                     "pmcoff_rec-wre_run-02", "pmcon_rec-wre_run-02"] 
                                for seq in ["t2tse", "t1tirm"]}
seq_to_relevant_subsequences["t2star"] = ["pmcoff_rec-wore_run-02", "pmcon_rec-wore_run-02"]
seq_to_relevant_subsequences["TRACEW"] = ["pmcoff_rec-wore_run-02", "pmcon_rec-wore_run-02"]



draw_boxplot_3(relevant_sequences_plot_3, seq_to_relevant_subsequences)