import numpy as np
import pandas as pd
from pathlib import Path
import pyvista as pv
from scipy.spatial.distance import euclidean

folder = Path('D:/VMWare Shared Folder/230915_all_conditions_mito_decon')

for file in folder.glob('*mitosurface.vtk'):
    print(file)
    mesh = pv.read(file)
    connected_components = mesh.connectivity()

    # Group cells by component label and calculate the sizes
    component_labels = connected_components['RegionId']
    unique_labels, label_counts = np.unique(
        component_labels, return_counts=True)
    component_sizes = dict(zip(unique_labels, label_counts))

    # Calculate the center of mass (centroid) for each connected component
    centroids = {}
    for label in unique_labels:
        mask = component_labels == label
        component_centers = mesh.extract_cells(mask)
        centroid = np.mean(component_centers.points, axis=0)
        centroids[label] = centroid

    # Convert the centroids to a PyVista Points object
    centroid_points = pv.PolyData(list(centroids.values()))

    # Find the center of mass for the entire image
    np_centroids = [np.array(centroid) for centroid in centroids.values()]
    cell_centroid = np.mean(np_centroids, axis=0)

    # Create a list to store component distance to CoM
    component_distance = []
    # Calculate distances from the image centroid to individual component centroids
    for centroid in np_centroids:
        distance = euclidean(cell_centroid, centroid)
        component_distance.append(distance)

    # Create a master list of filename, component, size and distance
    component_info = []
    for key, value in component_sizes.items():
        element_from_list = component_distance.pop(0)
        component_info.append([file.stem, key, value, element_from_list])

    df = pd.DataFrame(component_info, columns=[
                      'file', 'mito', 'size', 'distance'])
    df.to_csv('230915_all_conditions_decon/' +
              file.stem + '_mito_info.csv', index=False)
