# Evaluating the performance of markerless prospective motion correction and selective reacquisition in a general clinical protocol for brain MRI

Manuscript submitted to the Journal of Magnetic Resonance Imaging

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


## References:
[1] https://freesurfer.net/fswiki/recon-all, software version FreeSurfer v7.1.1 is used.

[2] Greve DN, Fischl B. Accurate and robust brain image alignment using boundary-based registration. NeuroImage. 2009;48(1):63-72. doi:10.1016/j.neuroimage.2009.06.060

[3] Reuter M, Rosas HD, Fischl B. Highly accurate inverse consistent registration: A robust approach. NeuroImage. 2010;53(4):1181-1196. doi:10.1016/j.neuroimage.2010.07.020
