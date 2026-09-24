'''
Reproduces a figure similar to the manuscript's Figure 3 ("Example images")
for a given subject: motion curves for the MPR (shake) and FLAIR (nod)
acquisitions, plus example image pairs (without PMC/reacquisition vs. with
PMC/reacquisition, where applicable) for T1_MPR, T2_FLAIR, T2_TSE, T1_TIRM,
T2* and the TRACE-weighted DWI.

No subject is named anywhere in the manuscript or repository for the
original Figure 3, so by default this reproduces it for two subjects
chosen out of the 7 subjects that have all 6 required sequences
(sub-02, 03, 07, 12, 13, 15, 17 -- limited by FLAIR/DWI, which were only
acquired in 10 of 22 subjects, and further by T2*, which is missing for
3 of those 10):
  - sub-02: the first of those 7 subjects by ID.
  - sub-03: the one with the largest mean Tenengrad improvement between
    "PMC off / reac off" and "PMC on / reac on", averaged over every
    sequence/run that has both variants (MPR nod, MPR shake, FLAIR nod,
    TSE nod, TIRM nod), computed from the metricsresults already
    produced by analysis_img_quality.py.

Usage:
    python3 make_figure3_example.py            # sub-02 and sub-03
    python3 make_figure3_example.py sub-07      # only sub-07
    python3 make_figure3_example.py sub-07 sub-12   # only those two
'''
import os
import sys
import datetime as dt

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

from motion_estimates import ExtractMotionParForScan

root = os.environ.get("MOCO_DATASET_PATH")
reg_dir = os.path.join(root, "derivatives/results/registrations")
out_dir = os.path.join(root, "derivatives/results/plots") + "/"
os.makedirs(out_dir, exist_ok=True)


def times_to_seconds(times):
    t0 = dt.datetime.strptime(times[0], "%H:%M:%S.%f")
    return np.array([(dt.datetime.strptime(t, "%H:%M:%S.%f") - t0).total_seconds()
                      for t in times])


def plot_motion(ax_t, ax_r, subj, name, seq_type, title):
    times, transl, rot = ExtractMotionParForScan(subj + '/', name, seq_type)
    t = times_to_seconds(times)

    for i, lab in enumerate(['R', 'A', 'S']):
        ax_t.plot(t, transl[:, i], lw=0.8, label=f'T_{lab}')
    ax_t.set_ylabel('Translation [mm]')
    ax_t.set_title(title)
    ax_t.legend(fontsize=7, loc='upper right')

    for i, lab in enumerate(['R', 'A', 'S']):
        ax_r.plot(t, rot[:, i], lw=0.8, label=f'R_{lab}')
    ax_r.set_ylabel('Rotation [deg]')
    ax_r.set_xlabel('Time from scan start [s]')
    ax_r.legend(fontsize=7, loc='upper right')


def mid_slice(img_path):
    img = nib.load(img_path).get_fdata()
    sh = np.shape(img)
    sl_dir = int(np.argmin(sh))
    indx = [slice(None)] * img.ndim
    indx[sl_dir] = int(sh[sl_dir] / 2)
    return np.squeeze(img[tuple(indx)])


def plot_image_pair(ax_off, ax_on, subj, tag, off_name, on_name, label, left_label, right_label):
    off_path = os.path.join(reg_dir, subj, tag, off_name)
    on_path = os.path.join(reg_dir, subj, tag, on_name)
    img_off = mid_slice(off_path)
    img_on = mid_slice(on_path)
    vmax = np.percentile(np.concatenate([img_off.ravel(), img_on.ravel()]), 99.5)

    ax_off.imshow(np.rot90(img_off), cmap='gray', vmin=0, vmax=vmax)
    ax_off.axis('off')
    ax_off.set_title(f'{label}\n{left_label}', fontsize=9, loc='left')

    ax_on.imshow(np.rot90(img_on), cmap='gray', vmin=0, vmax=vmax)
    ax_on.axis('off')
    ax_on.set_title(right_label, fontsize=9)


def make_figure(subj):
    fig = plt.figure(figsize=(14, 15.5))
    gs = fig.add_gridspec(7, 4, height_ratios=[1, 1, 2.2, 2.2, 2.2, 2.2, 2.2],
                           hspace=0.6, wspace=0.35, top=0.95, bottom=0.02)

    # (a) MPR shake motion curves, (b) FLAIR nod motion curves
    ax_a_t = fig.add_subplot(gs[0, 0:2])
    ax_a_r = fig.add_subplot(gs[1, 0:2])
    plot_motion(ax_a_t, ax_a_r, subj, 'mprage*pmcoff*_run-03_', 'run-03_mprage',
                'a: T1_MPR (SHAKE)')

    ax_b_t = fig.add_subplot(gs[0, 2:4])
    ax_b_r = fig.add_subplot(gs[1, 2:4])
    plot_motion(ax_b_t, ax_b_r, subj, 'flair*pmcoff*_run-02_', 'run-02_flair',
                'b: T2_FLAIR (NOD)')

    # (c)-(h) example image pairs
    subs = subj  # e.g. 'sub-02'

    panels = [
        ('c', 'T1_MPR', 'mprage',
         f'{subs}_acq-mpragepmcoff_rec-wore_run-03_T1w_moved.nii',
         f'{subs}_acq-mpragepmcon_rec-wre_run-03_T1w_moved.nii',
         'without PMC or reac', 'with PMC and reac'),
        ('d', 'T2_FLAIR', 'flair',
         f'{subs}_acq-flairpmcoff_rec-wore_run-02_FLAIR_moved.nii',
         f'{subs}_acq-flairpmcon_rec-wre_run-02_FLAIR_moved.nii',
         'without PMC or reac', 'with PMC and reac'),
        ('e', 'T2_TSE', 't2tse',
         f'{subs}_acq-t2tsepmcoff_rec-wore_run-02_T2w_moved.nii',
         f'{subs}_acq-t2tsepmcon_rec-wre_run-02_T2w_moved.nii',
         'without PMC or reac', 'with PMC and reac'),
        ('f', 'T1_TIRM', 't1tirm',
         f'{subs}_acq-t1tirmpmcoff_rec-wore_run-02_T1w_moved.nii',
         f'{subs}_acq-t1tirmpmcon_rec-wre_run-02_T1w_moved.nii',
         'without PMC or reac', 'with PMC and reac'),
        ('g', 'T2*', 't2star',
         f'{subs}_acq-t2starpmcoff_rec-wore_run-02_T2starw_moved.nii',
         f'{subs}_acq-t2starpmcon_rec-wore_run-02_T2starw_moved.nii',
         'without PMC', 'with PMC'),
        ('h', 'TRACEW DWI', 'TRACEWB1000',
         f'{subs}_acq-pmcoff_run-02_desc-TRACEWB1000_dwi_moved.nii',
         f'{subs}_acq-pmcon_run-02_desc-TRACEWB1000_dwi_moved.nii',
         'without PMC', 'with PMC'),
    ]

    grid_pos = [(2, 0, 2), (2, 2, 4), (3, 0, 2), (3, 2, 4), (4, 0, 2), (4, 2, 4)]
    # use rows 2-4 for c-h (2 panels per row, each panel = 2 columns split in half)
    for (label, seqname, tag, off_f, on_f, ltxt, rtxt), (row, c0, c1) in zip(panels, grid_pos):
        sub_gs = gs[row, c0:c1].subgridspec(1, 2, wspace=0.05)
        ax_off = fig.add_subplot(sub_gs[0, 0])
        ax_on = fig.add_subplot(sub_gs[0, 1])
        try:
            plot_image_pair(ax_off, ax_on, subj, tag, off_f, on_f,
                             f'{label}: {seqname}', ltxt, rtxt)
        except FileNotFoundError as e:
            ax_off.text(0.5, 0.5, f'missing:\n{e}', ha='center', va='center', fontsize=6, wrap=True)
            ax_off.axis('off')
            ax_on.axis('off')

    fig.suptitle(f'Figure-3-style example, {subj}', fontsize=16, y=0.995)

    out_path = out_dir + f'Fig3_Example_{subj}.png'
    fig.savefig(out_path, bbox_inches='tight', dpi=200)
    plt.close(fig)
    print('Saved', out_path)


if __name__ == "__main__":
    subjects = sys.argv[1:] if len(sys.argv) > 1 else ["sub-02", "sub-03"]
    for subj in subjects:
        make_figure(subj)
