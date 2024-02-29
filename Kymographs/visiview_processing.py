import numpy as np
from pathlib import Path
import tifffile

folder = input('Enter the folder path: ')
folder = Path(folder)
timestep = int(input('Enter the time interval between frames in seconds: '))

# generate a list of all the .nd files in the folder
filename_list = [str(pp) for pp in Path(folder).glob('*.nd')]
# file = open(filename + '.nd', 'r')

for filename in folder.glob('*.nd'):
    print(filename)
    file = open(filename)

    info_dict = {}
    for line in file:
        line = line.strip()
        if not line or line.startswith('#') or line == '"EndFile"':
            continue
        key, value = line.split(',', 1)
        key = key.strip('" ')
        value = value.strip('" ')
        if key in ['DoTimelapse', 'DoStage', 'DoWave', 'DoZSeries', 'WaveInFileName']:
            value = value == 'TRUE'

        info_dict[key] = value

    if bool(info_dict['DoTimelapse']):
        timepoints = int(info_dict['NTimePoints'])
    if bool(info_dict['DoStage']):
        stage_positions = int(info_dict['NStagePositions'])
    if bool(info_dict['DoWave']):
        wavelengths = int(info_dict['NWavelengths'])
    if bool(info_dict['DoZSeries']):
        z_steps = int(info_dict['NZSteps'])
        z_stepsize = float(info_dict['ZStepSize'])

    movie = np.zeros((timepoints, wavelengths, 2048, 2048), dtype=np.uint16)

    for s in range(1, stage_positions+1):
        for w in range(1, wavelengths+1):
            for t in range(1, timepoints+1):
                image_name = filename.stem + '_w' + \
                    str(w) + info_dict['WaveName' + str(w)
                                       ] + '_s{}_t{}.ome.tif'.format(s, t)
                image = tifffile.imread(folder / image_name)
                movie[t-1, w-1, :, :] = image
        tifffile.imwrite(
            folder.parent.joinpath(filename.stem + '_s' + str(s) + '.tif'),
            movie,
            imagej=True,
            resolution=(1/0.103, 1/0.103),
            metadata={'unit': 'um', 'axes': 'TCYX', 'finterval': timestep},
            dtype=movie.dtype,
            # shape=movie.shape,
        )
