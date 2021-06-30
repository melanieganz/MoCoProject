# RealNoiseMRI - MRI reconstruction challenge with realistic noise

[Challenge Website](https://realnoisemri.grand-challenge.org/Description/)  |  [MICCAI 2021 Satellite Event](https://miccai2021.org/en/MICCAI2021-CHALLENGES.html)

In recent years, there is a growing focus on the application of fast magnetic resonance imaging (MRI) based on prior knowledge. However, in most cases the generation of the under-sampled k-space data is artificial and not realistic, since patient motion corrupts the k-space differently than removing certain lines of k-space.

We organize a challenge, the RealNoiseMRI challenge, that aims at checking the robustness of fast MRI reconstruction methods to motion artifacts, by providing image and k-space datasets that consist of ground truth as well as motion degraded scans. We aim to closely mimic the ongoing challenge organized by Facebook at NeurIPS: https://fastmri.org/, but provide motion degraded data to realistically test robustness. 

More information about our challenge and sign up can be found on our [website](https://realnoisemri.grand-challenge.org/Description/).

## Available code:
* [recon_ismrmrd_dataset](https://github.com/ismrmrd/ismrmrd-python-tools/blob/master/recon_ismrmrd_dataset.py): useful script for loading h5 data as numpy array and performing a basic reconstruction using the Inverse Fourier Transform
* `Save_reconstruction.py`: helper functions for sorting the reconstructed image and saving it as .nii file
* `Evaluate_SSIM.py`: script for evaluating image quality metrics, i.e. SSIM, on the challenge test submissions
* `Rank_submissions.py`: script for ranking the test submissions based on the median rank profile for SSIM values
* `Utils.py`: utilitary functions
