# Evaluating the performance of markerless prospective motion correction and selective reacquisition in a general clinical protocol for brain MRI

Code for the manuscript titled "Evaluating the performance of markerless prospective motion correction and selective reacquisition in a general clinical protocol for brain MRI" (in submission). The presented code was adapted so that it can be run on the anonymized, published dataset ([ds004332](https://openneuro.org/datasets/ds004332) on OpenNeuro).
Furthermore, all metrics from the repository https://github.com/melanieganz/MoCoProject/tree/main/ImageQualityMetrics are incorporated into the scripts, as well, in order to enable re-calculation of all these metrics.

## Setup

1) **Get the dataset.** Download [ds004332](https://openneuro.org/datasets/ds004332). Every script below reads its input and writes its output relative to a single environment variable:
    ```bash
    export MOCO_DATASET_PATH=/path/to/ds004332-download/   # trailing slash required
    ```
    All pipeline outputs (registered images, metrics, plots, ...) are written under `$MOCO_DATASET_PATH/derivatives/results/`, alongside the dataset itself.

2) **Python environment.** `requirements.txt` is a pinned set of package versions verified to work on a plain Python 3.12 venv:
    ```bash
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    ```
    (`environment.yml` and `MoCoHealthy-spec-file.txt` are the original 2022 conda environment files; their pinned conda builds no longer resolve cleanly on current systems, so `requirements.txt` is the recommended path now.)

3) **FreeSurfer.** Needed by `analysis_img_quality.py` (`mri_vol2vol`, `bbregister`, `mri_robust_register`, `mri_binarize`) and `analysis_cort_thickness.py` (`recon-all`). The manuscript's own analysis used FreeSurfer v7.1.1 [1]; this code has also been verified against FreeSurfer 8.2.0. After installing, source its setup script (and make sure a license file is in place) before running anything:
    ```bash
    export FREESURFER_HOME=/path/to/freesurfer
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
    ```

4) **Headless plotting.** The plot-generating scripts call `plt.show()`; set `export MPLBACKEND=Agg` if running without a display.

## Code used for analysing the data:
* `analysis_img_quality.py`: registers each sequence to its ground-truth scan (via FreeSurfer's `bbregister`/`robust_register`) and calculates image quality metrics (SSIM, PSNR, Tenengrad, AES, gradient/image entropy, CoEnt) on the brain-masked images [1].
* `analysis_motion_data.py`: calculates motion metrics (RMS/median/maximum displacement) for each sequence and subject, statistically compares scan types, and optionally plots them.
* `analysis_cort_thickness.py`: generates thickness maps for each subject and each motion type / motion correction setting and fits general linear models to the thickness maps. **Not currently runnable against the OpenNeuro release**: it still depends on a longitudinal FreeSurfer processing scheme and directory layout from the original internal analysis pipeline that has no counterpart in the public dataset. This is a known, unresolved limitation, distinct from the path-portability fixes elsewhere in this repo.
* `img_quality_metrics.py`: functions for calculating the image quality metrics.
* `motion_estimates.py`: functions for loading the tracking data corresponding to a specific scan and for calculating motion metrics.
* `recon_register.py`: functions for running the FreeSurfer command `recon-all` and for registering images with `bbregister`[2] and `robust_register`[3].
* `statistical_tests.py`: functions to perform Wilcoxon signed rank tests.
* `utils.py`: plotting and statistics utility functions.
* `plot_generation_rewrite.py`: generates the manuscript's image-quality and ADC-histogram figures. This is the current, maintained plotting script.
* `generate_plots_manuscript.py`: the original plotting script `plot_generation_rewrite.py` was rewritten from. Superseded by it; kept for reference, not the recommended entry point.
* `make_figure3_example.py`: reproduces a figure similar to the manuscript's Figure 3 (motion curves plus example images with/without PMC and reacquisition) for a chosen subject. No such script existed in the original repository; see the script's docstring for how its default subjects were chosen.
* `raw_data_utils/`: scripts that need the *raw* (non-anonymized) scanner data rather than the OpenNeuro release -- see [raw_data_utils/README.md](raw_data_utils/README.md).

## Detailed instructions to rerun the analysis:
The analysis can be re-run in the following order. All steps read and write under `$MOCO_DATASET_PATH` (see Setup above).

1) Image quality assessment:
    Run the script `analysis_img_quality.py`. For a re-calculation of the metrics without redoing the registration -- the default, and the recommended way to run it against the OpenNeuro release, which already ships the registration transforms -- run with the following parameters:
    ```
    recon_all = False
    register = False
    apply_transform_bm = False
    apply_transform = True
    metrics = True
    show_bm_reg = False
    ```
    For redoing the FreeSurfer analysis and the registration from scratch, set `recon_all` and `register` to `True` as well (this reruns `recon-all` on every ground-truth scan and takes a long time). The script is resumable: on a rerun, subjects whose metrics already exist under `derivatives/results/metricsresults/` are skipped.

2) Analysis of subjects' motion:
    Run the script `analysis_motion_data.py` with `new_calc=True` to (re-)calculate the motion metrics on the anonymized dataset. Set `plot=True` as well to also generate the RMS/median/maximum-displacement boxplots (STILL, and NOD+SHAKE) as `Fig2_Motion_Boxplot_run-0{1,2}_*.png` under `derivatives/results/plots/`; if the metrics were already calculated in a previous run, set `new_calc=False` to skip straight to plotting.
    Note: `motion_estimates.py` reads scans' acquisition times from the BIDS JSON sidecars rather than the DICOM header.

3) Analysis of cortical thickness maps:
    Run the script `analysis_cort_thickness.py`. The motion data needs to be analysed first. **Not currently runnable against the OpenNeuro release** -- see the note above.

4) Generate the remaining plots for the manuscript:
    Run the script `plot_generation_rewrite.py`. It has no on/off flags: running it regenerates all four remaining figures (three image-quality boxplots and the ADC histogram) unconditionally, from whatever image-quality metrics and observer scores are present under `derivatives/`.

5) (Optional) Figure-3-style example images:
    Run `make_figure3_example.py [subject ...]` to reproduce a figure like the manuscript's Figure 3 (example images plus motion curves) for one or more subjects. Needs step 1's registered images to already exist.

Lastly, the script `raw_data_utils/get_scan_end.py` can be used to extract the scan end times again. This is optional, since we provide the extracted scan end times in `source/`. Unlike the rest of this pipeline, it needs the raw (non-anonymized) ismrmrd scanner data, which is available on request from PublicnEUro rather than from the OpenNeuro release -- see [raw_data_utils/README.md](raw_data_utils/README.md) for details.

## References:
[1] https://freesurfer.net/fswiki/recon-all, software version FreeSurfer v7.1.1 is used.

[2] Greve DN, Fischl B. Accurate and robust brain image alignment using boundary-based registration. NeuroImage. 2009;48(1):63-72. doi:10.1016/j.neuroimage.2009.06.060

[3] Reuter M, Rosas HD, Fischl B. Highly accurate inverse consistent registration: A robust approach. NeuroImage. 2010;53(4):1181-1196. doi:10.1016/j.neuroimage.2010.07.020
