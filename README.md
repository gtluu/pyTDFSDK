# pyTDFSDK

## About
This package is a Python wrapper for Bruker's TDF-SDK data access library to be used with other Python packages.

## Installation
This package can be installed to a Python virtual environment using `pip`. 
```
pip install git+https://github.com/gtluu/pyTDFSDK
```

## Example Usage
```python
import pandas as pd
from pyTDFSDK.init_tdf_sdk import init_tdf_sdk_api
from pyTDFSDK.classes import TsfData
from pyTDFSDK.tsf import tsf_read_line_spectrum_v2


# Initialize the TDF-SDK library using the packaged .dll or .so files for Windows or Linux, respectively.
dll = init_tdf_sdk_api()

data = TsfData(bruker_d_folder_name='some_timstof_data.d', tdf_sdk=dll)

# Get all spectra from a TSF file as a list of pd.DataFrames with columns for m/z and intensity in centroid mode. 
spectra_dfs = []
for index, row in data.analysis['Frames'].iterrows():
    mz_array, intensity_array = tsf_read_line_spectrum_v2(tdf_sdk=dll, handle=data.handle, frame_id=int(row['Id']))
    spectra_dfs.append(pd.DataFrame({'mz': mz_array, 'intensity': intensity_array}))

# Print number of spectra parsed from TSF file.
print(len(spectra_dfs))
# Print spectrum from first frame.
print(spectra_dfs[0])
```
