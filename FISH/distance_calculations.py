import numpy as np
from pathlib import Path
import pandas as pd
import bigfish
import bigfish.stack as stack
import bigfish.multistack as multistack


def resize_image(smaller_image, larger_image):
    rows_diff = larger_image.shape[0] - smaller_image.shape[0]
    cols_diff = larger_image.shape[1] - smaller_image.shape[1]

    new_arr = np.zeros(larger_image.shape, dtype=smaller_image.dtype)
    new_arr[:smaller_image.shape[0], :smaller_image.shape[1]] = smaller_image

    return new_arr


def select_last_slice(arr):
    '''For some reason napari decides to save a layer image as a 3D image instead of a 2D image. Therefore I have to remove the extra slices.
    I do thisby selecting the smallest dimension because for some reason it sometimes puts it first but also last.
    By comparing the first and last dimension I can determine which one is the smallest and then select the last slice.
    '''
    if arr.shape[0] < arr.shape[-1]:
        return arr[-1, :, :]
    else:
        return arr[:, :, -1]


def spot_distance_to_CoM(list_of_points, CoM, filename, dimensions, point_location):
    distance_to_CoM = []
    for point in list_of_points:
        distance = np.sqrt(
            (CoM[0] - point[dimensions-2]) ** 2 + (CoM[1] - point[dimensions-1]) ** 2)
        distance_to_CoM.append(
            (filename.stem, point[dimensions-2], point[dimensions-1], distance, point_location))
    return distance_to_CoM


folder = Path('data/231107_cells')

# load data
for file in folder.glob('*.tif'):
    print(file)
    mito_label = stack.read_image('output/' + file.stem + '_mito_label.tif')
    nuc_label = stack.read_image('output/' + file.stem + '_nuc_label.tif')
    cell_label = stack.read_image('output/' + file.stem + '_cell_label.tif')
    if len(cell_label.shape) > 2:
        cell_label = select_last_slice(cell_label)
    cell_label = cell_label.astype(np.int32)
    # match size of cell_label to size of the original image (via mito_label)
    cell_label = resize_image(cell_label, mito_label)
    spots = stack.read_array_from_csv(
        'output/' + file.stem + '_spots.csv', dtype=np.int64)
    random_spots = stack.read_array_from_csv(
        'output/' + file.stem + '_random_spots.csv', dtype=np.int64)

    # regular spot location identification
    spots_in_cell, spots_out_cell = multistack.identify_objects_in_region(
        cell_label, spots, ndim=3)
    spots_in_nuc, spots_out_nuc = multistack.identify_objects_in_region(
        nuc_label, spots_in_cell, ndim=3)
    spots_in_mito, spots_out_mito = multistack.identify_objects_in_region(
        mito_label, spots_out_nuc, ndim=3)
    # same for random spots
    rand_spots_in_cell, rand_spots_out_cell = multistack.identify_objects_in_region(
        cell_label, random_spots, ndim=2)
    rand_spots_in_nuc, rand_spots_out_nuc = multistack.identify_objects_in_region(
        nuc_label, rand_spots_in_cell, ndim=2)
    rand_spots_in_mito, rand_spots_out_mito = multistack.identify_objects_in_region(
        mito_label, rand_spots_out_nuc, ndim=2)

    # remove mito outside cell
    mito_label = mito_label*cell_label
    # calculate center of mass of mito (MTOC)
    center_of_mass = [np.average(indices)
                      for indices in np.where(mito_label == True)]

    # calculate distance of mito to the center of mass
    distances_mito = []
    for x, y in zip(*np.where(mito_label*cell_label)):
        distance = np.sqrt(
            (center_of_mass[0] - x) ** 2 + (center_of_mass[1] - y) ** 2)
        distances_mito.append((file.stem, distance))

    # calculate distance of mito FISH spots to center of mass
    distances_mito_FISH = spot_distance_to_CoM(
        spots_in_mito, center_of_mass, file, 3, 'mito')
    distances_mito_rand_FISH = spot_distance_to_CoM(
        rand_spots_in_mito, center_of_mass, file, 2, 'mito')

    # calculate distance of nucleus FISH spots to center of mass
    distances_nuc_FISH = spot_distance_to_CoM(
        spots_in_nuc, center_of_mass, file, 3, 'nuc')
    distances_nuc_rand_FISH = spot_distance_to_CoM(
        rand_spots_in_nuc, center_of_mass, file, 2, 'nuc')

    # calculate distance of all other FISH spots to center of mass
    distances_rest_FISH = spot_distance_to_CoM(
        spots_out_mito, center_of_mass, file, 3, 'out')
    distances_rest_rand_FISH = spot_distance_to_CoM(
        rand_spots_out_mito, center_of_mass, file, 2, 'out')

    spots_df = pd.DataFrame([[file.stem, center_of_mass[0], center_of_mass[1], 0, 'CoM']], columns=[
                            'file', 'x', 'y', 'distance', 'location'])
    spots_df = spots_df.append(pd.DataFrame(distances_mito_FISH, columns=[
                               'file', 'x', 'y', 'distance', 'location']))
    spots_df = spots_df.append(pd.DataFrame(distances_nuc_FISH, columns=[
                               'file', 'x', 'y', 'distance', 'location']))
    spots_df = spots_df.append(pd.DataFrame(distances_rest_FISH, columns=[
                               'file', 'x', 'y', 'distance', 'location']))

    rand_spots_df = pd.DataFrame([[file.stem, center_of_mass[0], center_of_mass[1], 0, 'CoM']], columns=[
        'file', 'x', 'y', 'distance', 'location'])
    rand_spots_df = rand_spots_df.append(pd.DataFrame(distances_mito_rand_FISH, columns=[
        'file', 'x', 'y', 'distance', 'location']))
    rand_spots_df = rand_spots_df.append(pd.DataFrame(distances_nuc_rand_FISH, columns=[
        'file', 'x', 'y', 'distance', 'location']))
    rand_spots_df = rand_spots_df.append(pd.DataFrame(distances_rest_rand_FISH, columns=[
        'file', 'x', 'y', 'distance', 'location']))

    mito_df = pd.DataFrame(distances_mito, columns=['file', 'distance'])

    pd.DataFrame.to_csv(spots_df, 'output/' +
                        file.stem + '_fish_spots_distances.csv')
    pd.DataFrame.to_csv(rand_spots_df, 'output/' +
                        file.stem + '_random_spots_distances.csv')
    pd.DataFrame.to_csv(mito_df, 'output/' + file.stem + '_mito_distances.csv')
