import os
import skimage.io
from skimage import data, color
from skimage.transform import rescale, resize, downscale_local_mean

def read_save_tif_lowres(path_morphology_focus, file_name, scale_method = ("rescale", "resize", "downscale_local_mean"),*args, **kwargs):
    fullres_img_tiff = skimage.io.imread(os.path.join(path_morphology_focus,file_name))
    if scale_method == "rescale":
        low_res = rescale(fullres_img_tiff, 0.1, anti_aliasing=True)
    elif scale_method == "resize":
        low_res = resize(fullres_img_tiff, (fullres_img_tiff.shape[0] // 10, fullres_img_tiff.shape[1] // 10), anti_aliasing=True)
    elif scale_method == "downscale_local_mean":
        low_res = downscale_local_mean(fullres_img_tiff, (4, 3))
    else:
        raise ValueError("Unknown scaling method: {}".format(scale_method))


    
    

    