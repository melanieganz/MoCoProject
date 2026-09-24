'''
Reproduces the manuscript's Figure 8 (motion-related cortical thickness
changes): FreeSurfer cross-sectional recon-all on the T1 MPR scans needed
per subject, resampling+smoothing cortical thickness to fsaverage, a
paired vertex-wise GLM (thickness ~ RMS motion) per condition against the
"Still without PMC" reference, FDR correction, and a composite figure in
the style of Fig. 8.

This corresponds to the paper's own description of that analysis exactly:
cortical thickness was generated with FreeSurfer's *cross-sectional*
recon-all stream (not the longitudinal base+long stream some earlier code
in this repo explored), using all T1 MPR scans, for all participants that
have them (all 22, by default -- T1 MPR is the one sequence acquired for
every subject).

Figure 8's five conditions, each compared against the same "Still without
PMC" reference:
    Without reacquisition: Still with PMC | Nod without PMC | Nod PMC
    With reacquisition:    Nod without PMC | Nod with PMC
(SHAKE is not part of this analysis; the paper's own text and Fig. 8 only
use STILL and NOD motion for the thickness/GLM analysis, even though other
code in this repo's history explored SHAKE and longitudinal streams too.)

Known simplification vs. the published method: the manuscript extracted the
"Still without PMC" reference's brain mask with FreeSurfer, then manually
corrected it, and reused that corrected mask via bbregister-space alignment
for the image-quality metrics (not for cortical thickness itself). This
script instead runs `recon-all -all` fully automatically for every scan,
including the reference, using FreeSurfer's own automated skull-stripping
throughout -- there is no manual QC/correction step. This only affects the
"Still without PMC" condition; every other condition needs its own
independent recon-all regardless, since each reflects that scan's own
motion artifacts.

Usage
-----
Only two things need to be supplied, both as environment variables (same
convention as the rest of this pipeline):

    export MOCO_DATASET_PATH=/path/to/ds004332-download/    # trailing slash
    export FREESURFER_HOME=/path/to/freesurfer

    python3 analysis_cort_thickness.py

Requires step 2 (analysis_motion_data.py, new_calc=True) to have already
been run, since the RMS motion values used as the GLM covariate are read
from its output under derivatives/results/Motion_Estimates*.

Optional flags:
    --subjects sub-01 sub-02 ...   only these subjects (default: all 22)
    --jobs N                       concurrent recon-all jobs (default 5)
    --threads-per-job N            OMP_NUM_THREADS per recon-all (default 2)
                                    (jobs * threads-per-job should not
                                    exceed the CPU budget you want to use;
                                    defaults to 5*2=10)
    --stage {recon,glm,plot,all}   run only part of the pipeline
                                    (default: all). recon-all is by far
                                    the slow part (hours per scan) and is
                                    resumable -- already-completed subjects
                                    are skipped on a rerun, so it's safe to
                                    interrupt and restart with --stage recon,
                                    then run --stage glm once it's done.

Output (all under $MOCO_DATASET_PATH/derivatives/results/):
    freesurfer_subjects/<sub>_<condition>/    recon-all output per scan
    cort_thickness/<condition>.fsgd           per-condition design file
    cort_thickness/{lh,rh}.<condition>.glmdir/  mri_glmfit output
    plots/Fig8_Cortical_Thickness.png         composite figure
'''
import argparse
import glob
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from utils import SortFiles  # noqa: E402


ALL_SUBJECTS = [f"sub-{i:02d}" for i in range(1, 23)]

# (condition_key, pmc, run, reac). reac is the BIDS 'rec-' tag: 'wore'
# (without reacquisition) or 'wre' (with). Still (run-01) only ever has
# 'wore', since reacquisition wasn't applied to scans without intentional
# motion.
REFERENCE = ("ref", "pmcoff", "run-01", "wore")

CONDITIONS = [
    ("still_pmc",    "pmcon",  "run-01", "wore"),
    ("nod_off",      "pmcoff", "run-02", "wore"),
    ("nod_on",       "pmcon",  "run-02", "wore"),
    ("nod_off_reac", "pmcoff", "run-02", "wre"),
    ("nod_on_reac",  "pmcon",  "run-02", "wre"),
]

CONDITION_TITLES = {
    "still_pmc":    "Still with PMC",
    "nod_off":      "Nod without PMC",
    "nod_on":       "Nod PMC",
    "nod_off_reac": "Nod without PMC\n(reac)",
    "nod_on_reac":  "Nod with PMC\n(reac)",
}


def nifti_path(root, sub, pmc, run, reac):
    return os.path.join(root, sub, 'anat',
                         f'{sub}_acq-mprage{pmc}_rec-{reac}_{run}_T1w.nii')


def fs_subject_id(sub, condition_key):
    return f'{sub}_{condition_key}'


# ---------------------------------------------------------------------------
# Motion metrics (RMS displacement), reused from analysis_motion_data.py's
# output. run-01 = STILL, run-02 = NOD (confirmed subject/run correspondence
# with the manuscript's own naming).
# ---------------------------------------------------------------------------

def load_rms_lookup(root):
    '''
    Returns {(sub, 'off'|'on', 'run-01'|'run-02'): RMS_displacement}, read
    from the mprage motion metrics analysis_motion_data.py already wrote
    under derivatives/results/Motion_Estimates*.
    '''
    lookup = {}
    for run in ('run-01', 'run-02'):
        pattern = os.path.join(root, 'derivatives/results',
                                f'Motion_EstimatesMotionMetrics_{run}', 'mprage_*.txt')
        files = glob.glob(pattern)
        if not files:
            raise FileNotFoundError(
                f"No motion metrics found matching {pattern}. Run "
                "analysis_motion_data.py (new_calc=True) first -- see "
                "step 2 of the top-level README.")
        f = SortFiles(files)[0]
        data = np.loadtxt(f, skiprows=1)
        if data.shape[0] != len(ALL_SUBJECTS):
            raise ValueError(f"{f} has {data.shape[0]} rows, expected {len(ALL_SUBJECTS)}")
        for i, sub in enumerate(ALL_SUBJECTS):
            lookup[(sub, 'off', run)] = data[i, 0]   # RMS_Off
            lookup[(sub, 'on', run)] = data[i, 3]    # RMS_On
    return lookup


def rms_for_condition(rms_lookup, sub, pmc, run):
    status = 'off' if pmc == 'pmcoff' else 'on'
    # nod_off_reac/nod_on_reac share the same physical motion (and hence
    # the same tracked RMS) as nod_off/nod_on -- reacquisition changes the
    # reconstruction, not the motion that was tracked.
    return rms_lookup[(sub, status, run)]


# ---------------------------------------------------------------------------
# FreeSurfer command execution
# ---------------------------------------------------------------------------

def run_fs(cmd, freesurfer_home, subjects_dir, omp_threads=1, log_file=None):
    '''Runs a FreeSurfer command with the environment SetUpFreeSurfer.sh would set.'''
    setup = (f'export FREESURFER_HOME={freesurfer_home}; '
              f'source {freesurfer_home}/SetUpFreeSurfer.sh > /dev/null 2>&1; '
              f'export SUBJECTS_DIR={subjects_dir}; '
              f'export OMP_NUM_THREADS={omp_threads}; ')
    full_cmd = ['bash', '-c', setup + cmd]
    if log_file:
        with open(log_file, 'a') as f:
            f.write(f'\n=== {cmd} ===\n')
            f.flush()
            return subprocess.run(full_cmd, stdout=f, stderr=subprocess.STDOUT).returncode
    return subprocess.run(full_cmd).returncode


def ensure_fsaverage(freesurfer_home, subjects_dir):
    target = os.path.join(subjects_dir, 'fsaverage')
    if not os.path.exists(target):
        src = os.path.join(freesurfer_home, 'subjects', 'fsaverage')
        os.makedirs(subjects_dir, exist_ok=True)
        os.symlink(src, target)


# ---------------------------------------------------------------------------
# Stage 1: recon-all (the slow part)
# ---------------------------------------------------------------------------

def recon_all_job(root, freesurfer_home, subjects_dir, sub, condition_key,
                   pmc, run, reac, threads, log_dir):
    fs_id = fs_subject_id(sub, condition_key)
    done_marker = os.path.join(subjects_dir, fs_id, 'scripts', 'recon-all.done')
    if os.path.exists(done_marker):
        return (fs_id, 'skipped (already done)')

    nifti = nifti_path(root, sub, pmc, run, reac)
    if not os.path.exists(nifti):
        return (fs_id, f'ERROR: missing input {nifti}')

    log_file = os.path.join(log_dir, f'{fs_id}.log')
    cmd = f'recon-all -i {nifti} -s {fs_id} -all -parallel'
    rc = run_fs(cmd, freesurfer_home, subjects_dir, omp_threads=threads, log_file=log_file)
    ok = os.path.exists(done_marker)
    if ok:
        return (fs_id, 'done')
    return (fs_id, f'FAILED (exit {rc}, see {log_file})')


def stage_recon(root, freesurfer_home, subjects_dir, subjects, jobs, threads_per_job):
    log_dir = os.path.join(subjects_dir, '..', 'recon_all_logs')
    log_dir = os.path.normpath(log_dir)
    os.makedirs(log_dir, exist_ok=True)
    ensure_fsaverage(freesurfer_home, subjects_dir)

    tasks = []
    for sub in subjects:
        for key, pmc, run, reac in [REFERENCE] + CONDITIONS:
            tasks.append((sub, key, pmc, run, reac))

    total_cores = jobs * threads_per_job
    print(f'{len(tasks)} recon-all jobs queued '
          f'({jobs} concurrent x {threads_per_job} threads/job = {total_cores} cores)')

    results = []
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futures = {
            ex.submit(recon_all_job, root, freesurfer_home, subjects_dir,
                      sub, key, pmc, run, reac, threads_per_job, log_dir): (sub, key)
            for sub, key, pmc, run, reac in tasks
        }
        for fut in as_completed(futures):
            fs_id, status = fut.result()
            print(f'[{time.strftime("%H:%M:%S")}] {fs_id}: {status}', flush=True)
            results.append((fs_id, status))

    failed = [r for r in results if 'FAILED' in r[1] or 'ERROR' in r[1]]
    print(f'\nrecon-all stage done: {len(results)-len(failed)}/{len(results)} ok.')
    if failed:
        print(f'{len(failed)} FAILED/missing -- see logs under {log_dir}:')
        for fs_id, status in failed:
            print(f'  {fs_id}: {status}')
    return len(failed) == 0


# ---------------------------------------------------------------------------
# Stage 2: group GLM analysis
# ---------------------------------------------------------------------------

def write_fsgd(condition_key, pmc, run, subjects, rms_lookup, out_dir):
    lines = []
    for sub in subjects:
        rms_ref = rms_for_condition(rms_lookup, sub, REFERENCE[1], REFERENCE[2])
        rms_cond = rms_for_condition(rms_lookup, sub, pmc, run)
        lines.append(f'Input {fs_subject_id(sub, "ref")} Main {rms_ref}')
        lines.append(f'Input {fs_subject_id(sub, condition_key)} Main {rms_cond}')

    header = 'GroupDescriptorFile 1\nMeasurementName thickness\nClass Main\nVariables RMS'
    fsgd_path = os.path.join(out_dir, f'{condition_key}.fsgd')
    with open(fsgd_path, 'w') as f:
        f.write(header + '\n')
        f.write('\n'.join(lines) + '\n')
    return fsgd_path


def write_contrast(out_dir):
    '''
    A single-class FSGD with one continuous variable (RMS) has a 2-column
    design matrix: [offset, RMS-slope]. This contrast tests the RMS slope.
    '''
    path = os.path.join(out_dir, 'Avg-thickness-RMS-Cor.mtx')
    with open(path, 'w') as f:
        f.write('0 1\n')
    return path


def stage_glm(root, freesurfer_home, subjects_dir, subjects, rms_lookup, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    ensure_fsaverage(freesurfer_home, subjects_dir)
    contrast_path = write_contrast(out_dir)

    sig_fdr_paths = {}
    for condition_key, pmc, run, reac in CONDITIONS:
        print(f'--- {condition_key} ---', flush=True)
        fsgd_path = write_fsgd(condition_key, pmc, run, subjects, rms_lookup, out_dir)

        hemi_fdr = {}
        for hemi in ('lh', 'rh'):
            thick00 = os.path.join(out_dir, f'{hemi}.{condition_key}.thickness.00.mgh')
            thick10 = os.path.join(out_dir, f'{hemi}.{condition_key}.thickness.10.mgh')
            glmdir = os.path.join(out_dir, f'{hemi}.{condition_key}.glmdir')

            run_fs(f'mris_preproc --fsgd {fsgd_path} --target fsaverage --hemi {hemi} '
                   f'--meas thickness --out {thick00}', freesurfer_home, subjects_dir)
            run_fs(f'mri_surf2surf --hemi {hemi} --s fsaverage --sval {thick00} '
                   f'--fwhm 10 --cortex --tval {thick10}', freesurfer_home, subjects_dir)
            run_fs(f'mri_glmfit --y {thick10} --fsgd {fsgd_path} --C {contrast_path} '
                   f'--surf fsaverage {hemi} --cortex --glmdir {glmdir}',
                   freesurfer_home, subjects_dir)

            sig = os.path.join(glmdir, 'Avg-thickness-RMS-Cor', 'sig.mgh')
            sig_fdr = os.path.join(glmdir, 'Avg-thickness-RMS-Cor', 'sig_fdr.mgh')
            hemi_fdr[hemi] = (sig, sig_fdr)

        # FDR correction across both hemispheres together, as in the manuscript
        (lh_sig, lh_fdr), (rh_sig, rh_fdr) = hemi_fdr['lh'], hemi_fdr['rh']
        run_fs(f'mri_fdr --i {lh_sig} nomask {lh_fdr} --i {rh_sig} nomask {rh_fdr} --fdr 0.05',
               freesurfer_home, subjects_dir)

        sig_fdr_paths[condition_key] = {'lh': lh_fdr, 'rh': rh_fdr}

    return sig_fdr_paths


# ---------------------------------------------------------------------------
# Stage 3: composite figure (style of Fig. 8a)
# ---------------------------------------------------------------------------

def stage_plot(subjects_dir, out_dir, out_path):
    try:
        from nilearn import plotting as nlp
        import matplotlib.pyplot as plt
    except ImportError:
        print('nilearn is not installed -- skipping figure generation.\n'
              'Install it with: pip install nilearn\n'
              'The underlying sig_fdr.mgh significance maps are still available '
              f'under {out_dir}/{{lh,rh}}.<condition>.glmdir/Avg-thickness-RMS-Cor/, '
              'and can be viewed with e.g. freeview.')
        return

    fsaverage_surf = os.path.join(subjects_dir, 'fsaverage', 'surf')
    n_cond = len(CONDITIONS)
    fig, axes = plt.subplots(2, n_cond, figsize=(3.2 * n_cond, 6.4),
                              subplot_kw={'projection': '3d'})

    for col, (condition_key, *_rest) in enumerate(CONDITIONS):
        sig_fdr = os.path.join(out_dir, f'lh.{condition_key}.glmdir',
                                'Avg-thickness-RMS-Cor', 'sig_fdr.mgh')
        mesh = os.path.join(fsaverage_surf, 'lh.inflated')
        bg = os.path.join(fsaverage_surf, 'lh.curv')

        for row, view in enumerate(['lateral', 'medial']):
            ax = axes[row, col]
            if os.path.exists(sig_fdr):
                nlp.plot_surf_stat_map(
                    mesh, sig_fdr, hemi='left', view=view, bg_map=bg,
                    threshold=1.301,  # -log10(0.05)
                    cmap='cold_hot', colorbar=(col == n_cond - 1 and row == 0),
                    axes=ax, figure=fig,
                )
            else:
                ax.text2D(0.5, 0.5, 'missing', ha='center', transform=ax.transAxes)
                ax.axis('off')
            if row == 0:
                ax.set_title(CONDITION_TITLES[condition_key], fontsize=10)

    fig.suptitle('Cortical thickness vs. RMS motion (left hemisphere, FDR<0.05)',
                  fontsize=13)
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print('Saved', out_path)


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--subjects', nargs='+', default=ALL_SUBJECTS,
                         help='subjects to include (default: all 22)')
    parser.add_argument('--jobs', type=int, default=5,
                         help='concurrent recon-all jobs (default 5)')
    parser.add_argument('--threads-per-job', type=int, default=2,
                         help='OMP_NUM_THREADS per recon-all job (default 2; '
                              'jobs*threads-per-job defaults to 10 cores total)')
    parser.add_argument('--stage', choices=['recon', 'glm', 'plot', 'all'], default='all',
                         help='run only part of the pipeline (default: all)')
    args = parser.parse_args()

    root = os.environ.get('MOCO_DATASET_PATH')
    freesurfer_home = os.environ.get('FREESURFER_HOME')
    if not root:
        sys.exit('MOCO_DATASET_PATH is not set.')
    if not freesurfer_home:
        sys.exit('FREESURFER_HOME is not set.')

    subjects_dir = os.path.join(root, 'derivatives/results/freesurfer_subjects/')
    glm_dir = os.path.join(root, 'derivatives/results/cort_thickness/')
    plot_dir = os.path.join(root, 'derivatives/results/plots/')
    os.makedirs(subjects_dir, exist_ok=True)
    os.makedirs(glm_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    if args.stage in ('recon', 'all'):
        ok = stage_recon(root, freesurfer_home, subjects_dir, args.subjects,
                          args.jobs, args.threads_per_job)
        if not ok and args.stage == 'all':
            sys.exit('Some recon-all jobs failed -- fix and rerun with --stage recon '
                      'before continuing to --stage glm.')

    if args.stage in ('glm', 'all'):
        rms_lookup = load_rms_lookup(root)
        stage_glm(root, freesurfer_home, subjects_dir, args.subjects, rms_lookup, glm_dir)

    if args.stage in ('plot', 'all'):
        stage_plot(subjects_dir, glm_dir, os.path.join(plot_dir, 'Fig8_Cortical_Thickness.png'))


if __name__ == '__main__':
    main()
