import tifffile
import matplotlib.pyplot as plt
import os
import seaborn as sns
from skimage.transform import resize
import argparse
import numpy as np
import skimage

def raed_write_xinum_image(file_name,path_morphology_focus,scale,ch,prefix):
    file_path = os.path.join(path_morphology_focus, file_name)
    fullres_multich_img = tifffile.imread(file_path, is_ome=True, level=0, aszarr=False)
    fullres_multich_img_ch = fullres_multich_img[ch,:,:]
    image_resized = resize(fullres_multich_img_ch, (fullres_multich_img_ch.shape[0] // scale, fullres_multich_img_ch.shape[1] // scale), anti_aliasing=True)
    image_scaled = ((image_resized /image_resized.max())  * 255).astype(np.uint8)
    skimage.io.imsave(os.path.join(path_morphology_focus, f'_mean_{prefix}' + '.png'), image_scaled)

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description='sbatch create')
    parser.add_argument('--file_name', dest='file_name', type=str, required=True)
    parser.add_argument('--path_morphology_focus', dest='path_morphology_focus', type=str, required=True)
    parser.add_argument('--scale', dest='scale', type=int, required=True)
    parser.add_argument('--ch', dest='ch', type=int, required=True)
    parser.add_argument('--prefix', dest='prefix', type=str, required=True)
    args = parser.parse_args()
    raed_write_xinum_image(file_name = args.file_name,path_morphology_focus = args.path_morphology_focus,scale = args.scale, ch = args.ch, prefix = args.prefix)
