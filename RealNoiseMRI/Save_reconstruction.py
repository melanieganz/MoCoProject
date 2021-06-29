'''
Script for sorting the reconstructed image and saving it as .nii file
'''

import numpy as np
import nibabel as nib
from ismrmrdtools import show


def ReOrderImg(img, show_img=True):
    '''
    Re-order the slices, since they were acquired in an interleaved mode.

    Parameters
    ----------
    img : np.array
        Image reconstructed from raw data. Expects raw data to be loaded and
        image to be reconstructed as performed in: 
        https://github.com/ismrmrd/ismrmrd-python-tools/blob/master/recon_ismrmrd_dataset.py.
        In especially the image is expected to have the following dimensions:
        (nreps, ncontrasts, nslices, ncoils, eNz, eNy, rNx).
    show_img : bool, optional
        Whether the image and sorted image should be shown as overview. The 
        default is True.

    Returns
    -------
    img_res : np.array
        sorted image with the dimensions (eNy, rNx, nslices). Axes 0 and 1 are
        flipped to ensure correct right-left and anterior-posterior sorting.

    '''
    img = img[0,0,:,0,:,:]
    num = np.shape(img)[0]
    if (num%2) == 0:
        # even number of slices in img:
        indices = []
        for i in range(int(num/2)):
            indices.append(i)
            indices.append(i+int(num/2))
    else:
        # odd number of slices in img:
        n = int(num/2)
        indices = []
        for i in range(n):
            indices.append(i)
            indices.append(i+n+1)
        indices.append(n)
    
    img_sorted = img[indices]
    
    if show_img:
        # show the image first unsorted and then sorted
        show.imshow(np.squeeze(img))
        show.imshow(np.squeeze(img_sorted))
    
    # re-sort the image from axes (0,1,2) into (1,2,0):
    sh = np.shape(img_sorted)
    img_res = np.zeros((sh[1], sh[2], sh[0]))
    for i in range(sh[0]):
        img_res[:,:,i] = img_sorted[i]
    
    # flip axes 0 and 1, so that right-left and anterior-posterior are sorted 
    # correctly:
    img_res = img_res[::-1,::-1]
    
    return img_res    

 
def Save_Nifti(img, vox2ras, filename):
    '''
    Save the image as nifti file. Use vox2ras transform as affine transform 
    to RAS space for the nifti.

    Parameters
    ----------
    img : np.array
        Image to be saved.
    vox2ras : np.array
        Affine transform from voxel space to RAS space, which needs to be
        specified for Nifti files.
    filename : str
        Filename under which the image should be saved (without -nii extension).

    Returns
    -------
    None.

    '''
    
    nii = nib.Nifti1Image(img, vox2ras)
    nib.save(nii, filename+'.nii')



# test the functions for an exemplary image:
img = np.load('Example_03_t1.npy')
vox2ras = np.loadtxt('vox2ras_nod_t1_off.txt')
filename = 'Example_03_t1'
    
    
img_sorted = ReOrderImg(img)
Save_Nifti(img_sorted, vox2ras, filename)    
    
    
