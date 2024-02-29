import numpy as np
import random
from scipy.ndimage import binary_fill_holes
from tifffile import imread
from skimage import filters
from pathlib import Path
import bigfish
import bigfish.stack as stack
import bigfish.detection as detection


def select_last_slice(arr):
    '''For some reason napari decides to save a layer image as a 3D image instead of a 2D image. Therefore I have to remove the extra slices.
    I do thisby selecting the smallest dimension because for some reason it sometimes puts it first but also last.
    By comparing the first and last dimension I can determine which one is the smallest and then select the last slice.
    '''
    if arr.shape[0] < arr.shape[-1]:
        return arr[-1, :, :]
    else:
        return arr[:, :, -1]


folder = Path('data/231107_cells')
output = 'output_single_slice/'
voxel_size = (110, 110)  # in nanometer (one value per dimension zyx)
spot_size = (150, 150)

# load data
for file in folder.glob('*.tif'):
    print(file)
    image = imread(file)
    rna = image[:, 0, :, :]
    mito = image[:, 1, :, :]
    dapi = image[:, 2, :, :]

    # find most in focus mito slice
    mito_focus = stack.compute_focus(mito)
    mito_focus_index = stack.get_in_focus_indices(mito_focus, 1)

    mito_2d = mito[mito_focus_index[0]]
    dapi_2d = dapi[mito_focus_index[0]]
    rna_2d = rna[mito_focus_index[0]]

    # segment mitochondria
    mito_gauss = filters.gaussian(mito_2d, sigma=2)
    mito_label = mito_gauss > filters.threshold_li(mito_gauss)
    mito_label = np.multiply(mito_label, 1)

    # detect RNA spots and do subpixel fitting
    spots, threshold = detection.detect_spots(
        images=rna_2d,
        # threshold=150,  # only for cells from 231107!!
        return_threshold=True,
        voxel_size=voxel_size,
        spot_radius=spot_size)

    spots_subpixel_rna = detection.fit_subpixel(
        image=rna_2d,
        spots=spots,
        voxel_size=voxel_size,
        spot_radius=spot_size)

    cell_label = stack.read_image(
        'output/' + file.stem + '_cell_label.tif')
    if len(cell_label.shape) > 2:
        cell_label = select_last_slice(cell_label)

    # random points as control
    # Set the size of the image (x and y dimensions)
    x = cell_label.shape[0]  # Replace with your desired width
    y = cell_label.shape[1]  # Replace with your desired height
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
    stack.save_data_to_csv(
        spots_subpixel_rna, output + file.stem + '_spots.csv')
    stack.save_data_to_csv(random_points, output +
                           file.stem + '_random_spots.csv')
