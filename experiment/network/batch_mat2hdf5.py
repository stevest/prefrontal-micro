import numpy as np
#import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import h5py

print(__name__)
print(__name__ is "__main__")

if __name__ == "__main__":
    #===================================================================================================================
    # Have here the convertion from mat files to hdf5 files:
    print('Initializing script')
    root_dir = Path(r'/home/cluster/stefanos/Documents/Glia')
    #root_dir = Path(r'C:\Users\stefanos\Documents\analysis2\analysis\test')
    print('Locating files to be transformed...')
    files = list(root_dir.glob('**/vsoma.mat'))
    nfiles = len(files)
    print(f'{nfiles} files were found!')
    #fig, ax = plt.subplots()
    for i, file in enumerate(files):
        print(f'{i}/{nfiles}: Transforming file {file}')
        with h5py.File(file, 'r') as f:
            #list(f.keys())
            data = f['vsoma']
            #data = dset.value
            ncells = data.shape[1]

            ref = data[0, 0]
            deref = f[ref]
            rawdata = np.array(deref)
            nsamples = rawdata.size

            vsoma = np.empty([ncells, nsamples], dtype=float)

            for cell in range(ncells):
                ref = data[0, cell]
                deref = f[ref]
                rawdata = np.array(deref)
                #print(f'{cell} size: {rawdata.shape}')
                vsoma[cell][:] = rawdata
                #ax.plot(rawdata[0, ::10])
                #plt.pause(0.5)

            filename = str(file).replace('.mat', '.hdf5')
            df = pd.DataFrame(vsoma)
            df.to_hdf(filename, key='vsoma', mode='w')

    print('success!')
