import numpy as np
import random
from scipy.ndimage import binary_fill_holes
from tifffile import imread
from skimage import filters
from pathlib import Path
import bigfish
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.segmentation as segmentation

folder = Path('data/231107_cells')
output = 'output/'
voxel_size = (200, 110, 110)  # in nanometer (one value per dimension zyx)
spot_size = (250, 150, 150)

# load pretrained model for nucleus segmentation
model_nuc = segmentation.unet_3_classes_nuc()

# load data
for file in folder.glob('*.tif'):
    print(file)
    image = imread(file)
    rna = image[:, 0, :, :]
    mito = image[:, 1, :, :]
    dapi = image[:, 2, :, :]

    # max intensity projections
    mito_mip = stack.maximum_projection(mito)
    dapi_mip = stack.maximum_projection(dapi)
    rna_mip = stack.maximum_projection(rna)

    # segment mitochondria
    mito_gauss = filters.gaussian(mito_mip, sigma=2)
    mito_label = mito_gauss > filters.threshold_li(mito_gauss)
    mito_label = np.multiply(mito_label, 1)

    # segment the nucleus
    nuc_label = segmentation.apply_unet_3_classes(
        model_nuc, dapi_mip, target_size=256, test_time_augmentation=True)
    # fill holes in the nucleus segmentation
    nuc_label = binary_fill_holes(nuc_label)
    nuc_label = np.multiply(nuc_label, 1)

    # detect RNA spots and do subpixel fitting
    spots, threshold = detection.detect_spots(
        images=rna,
        # threshold=150, #only for cells from 231107!!
        return_threshold=True,
        voxel_size=voxel_size,
        spot_radius=spot_size)

    spots_subpixel_rna = detection.fit_subpixel(
        image=rna,
        spots=spots,
        voxel_size=voxel_size,
        spot_radius=spot_size)

    # random points as control
    # Set the size of the image (x and y dimensions)
    x = mito_mip.shape[0]  # Replace with your desired width
    y = mito_mip.shape[1]  # Replace with your desired height
    # Set the number of random points you want to generate
    num_points = spots.shape[0]  # Replace with your desired number of points
    # Create a list to store the random points
    random_points = []
    # Generate random points and add them to the list
    for _ in range(num_points):
        # Generate random x and y coordinates within the image dimensions
        random_x = random.randint(0, x - 1)
        random_y = random.randint(0, y - 1)

        # Append the point as a tuple (x, y) to the list
        random_points.append((random_x, random_y))

    random_points = np.array(random_points)

    # save segmentations and spots
    stack.save_image(mito_label, output + file.stem + '_mito_label.tif')
    stack.save_image(nuc_label, output + file.stem + '_nuc_label.tif')
    stack.save_data_to_csv(
        spots_subpixel_rna, output + file.stem + '_spots.csv')
    stack.save_data_to_csv(random_points, output +
                           file.stem + '_random_spots.csv')
