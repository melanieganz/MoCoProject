# Evaluating the performance of markerless prospective motion correction and selective reacquisition in a general clinical protocol for brain MRI

Code for the manuscript soon available on neurolivre.org (in preparation). The presented code was adapted so that it can be run on the anonymized, published dataset.
Furthermore, all metrics from the repository https://github.com/melanieganz/MoCoProject/tree/main/ImageQualityMetrics are incorporated into the scripts, as well, in order to enable re-calculation of all these metrics.

## Code used for analysing the data:
* `analysis_cort_thickness.py`: generates thickness maps for each subject and each motion tpye / motion correction setting and fits general linear models to the thickness maps
* `analysis_img_quality.py`: runs the groundtruth MPRAGE through FreeSurfers recon-all stream [1], performs image registration and calculates image quality metrics on the brain masked images
* `analysis_motion_data.py`: calculates motion metrics for each sequence and each subject and statistically compares different scan types
* `img_quality_metrics.py`: functions for calculating the image quality metrics
* `motion_estimates.py`: functions for loading the tracking data corresponding to a specific scan and for calculating motion metrics
* `recon_register.py`: functions for running the FreeSurfer command `recon-all` and for registering images with `bbregister`[2] and `robust_register`[3]
* `statistical_tests.py`: functions to perform Wilcoxon signed rank tests
* `utils.py`: utility functions
* `generate_plots_manuscript.py`: script to generate all plots of the manuscript

## Detailed instructions to rerun the analysis:
The analysis can be re-run in the following order:
1) Image quality assessment:
    Run the function Full_Analyis() in the script 'analysis_img_quality.py'. For a re-calculation of the metrics without redoing the registration please run with the following parameters:
    recon_all = False
    register = False
    apply_transform_bm = False
    apply_transform = True
    metrics = True
    show_bm_reg = False 
    For redoing the FreeSurfer analyis and the registration, all parameters need to be set to 'True'.

2) Analysis of subject's motion:
    Run the script 'analysis_motion_data.py' with the options new_calc=True and plot=False for re-calculating the motion metrics on the anonymized data set. The option plot=True does not need to be run, since the output plot will be generated in step 4.
    Note: The script 'motion_estimates.py' has been changed so that the scans' acquisition times ar not read from the DICOM header but from the JSON file.
    
3) Analyis of cortical thickness maps:
    Run the script 'analysis_cort_thickness.py'. The motion data needs to be analysed first.
    
4) Generate plots for manuscript:
    Run the script 'generate_plots_manuscript.py' with the following parameters:
    plot_motion = True
    plot_still = True
    plot_nod = True
    plot_DWI = True
    calc_ADC_hist = True or False (depending on whether the metrics for the ADC histograms should be re-calculated).
    
Lastly, the script 'get_scan_end.py' can be used to extract the scan end times again. This is optional, since we provide the extracted scan end times in ../TCLData.

## References:
[1] https://freesurfer.net/fswiki/recon-all, software version FreeSurfer v7.1.1 is used.

[2] Greve DN, Fischl B. Accurate and robust brain image alignment using boundary-based registration. NeuroImage. 2009;48(1):63-72. doi:10.1016/j.neuroimage.2009.06.060

[3] Reuter M, Rosas HD, Fischl B. Highly accurate inverse consistent registration: A robust approach. NeuroImage. 2010;53(4):1181-1196. doi:10.1016/j.neuroimage.2010.07.020
