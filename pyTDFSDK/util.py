import numpy as np
from ctypes import POINTER, c_double
from pyTDFSDK.error import *


def call_conversion_func(tdf_sdk, handle, frame_id, input_data, func):
    """
    Decorator for other conversion related functions.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_open() or
        pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param input_data: Array to be converted.
    :type input_data: numpy.array
    :param func: Conversion function from TDF-SDK to pass through decorator.
    :type func: function
    :return: Array of converted values equal in size to the input_data array.
    :rtype: numpy.array
    """
    if type(input_data) is np.ndarray and input_data.dtype == np.float64:
        # already "native" format understood by DLL -> avoid extra copy
        in_array = input_data
    else:
        # convert data to format understood by DLL
        in_array = np.array(input_data, dtype=np.float64)
    cnt = len(in_array)
    out = np.empty(shape=cnt, dtype=np.float64)
    success = func(handle,
                   frame_id,
                   in_array.ctypes.data_as(POINTER(c_double)),
                   out.ctypes.data_as(POINTER(c_double)),
                   cnt)
    if success == 0:
        throw_last_timsdata_error(tdf_sdk)
    return out


def get_encoding_dtype(encoding):
    """
    Use "encoding" command line parameter to determine numpy dtype.

    :param encoding: Encoding command line parameter, either "64" or "32".
    :type encoding: int
    :return: Numpy dtype, either float64 or float32
    :rtype: numpy.dtype
    """
    if encoding == 32:
        return np.float32
    elif encoding == 64:
        return np.float64


def get_centroid_status(mode, exclude_mobility=None):
    """
    Use "mode" command line parameter to determine whether output data is centroided in psims compatible format.

    :param mode: Mode command line parameter, either "profile", "centroid", or "raw".
    :type mode: str
    :param exclude_mobility: Whether to include mobility data in the output files, defaults to None.
    :type exclude_mobility: bool | None
    :return: Dictionary containing standard spectrum data.
    :return: Tuple of (centroided status (bool), exclude_mobility status (bool))
    :rtype: tuple[bool]
    """
    if mode == 'profile':
        centroided = False
        exclude_mobility = True
    elif mode == 'centroid' or mode == 'raw':
        centroided = True
    return centroided, exclude_mobility


def get_maldi_coords(data, maldiframeinfo_dict):
    """
    Get tuple of MALDI coordinates from analysis.tsf/analysis.tdf metadata.

    :param data: Object containing metadata from analysis.tsf/analysis.tdf database.
    :type data: pyTDFSDK.classes.TsfData | pyTDFSDK.classes.TdfData
    :param maldiframeinfo_dict: A row from the MaldiFrameInfo table in analysis.tsf/analysis.tdf database.
    :type maldiframeinfo_dict: dict
    :return: x-y (or x-y-z if available) coordinates for the current spectrum.
    :rtype: tuple[int]
    """
    if data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'SingleSpectra':
        coords = maldiframeinfo_dict['SpotName']
    elif data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'Imaging':
        coords = [int(maldiframeinfo_dict['XIndexPos']), int(maldiframeinfo_dict['YIndexPos'])]
        if 'ZIndexPos' in data.analysis['MaldiFrameInfo'].columns:
            coords.append(int(maldiframeinfo_dict['ZIndexPos']))
        coords = tuple(coords)
    return coords


def bin_profile_spectrum(mz_array, intensity_array, profile_bins, mz_encoding):
    """
    Bin profile mode spectrum into N number of bins.

    :param mz_array: Array containing m/z values.
    :type mz_array: numpy.array
    :param intensity_array: Array containing intensity values.
    :type intensity_array: numpy.array
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param mz_encoding: m/z encoding command line parameter, either "64" or "32".
    :type mz_encoding: int
    :return: Tuple of binned_mz_array (np.array) and binned_intensity_array (np.array).
    :rtype: tuple[numpy.array]
    """
    mz_acq_range_lower = float(mz_array[0])
    mz_acq_range_upper = float(mz_array[-1])
    bins = np.linspace(mz_acq_range_lower, mz_acq_range_upper, profile_bins, dtype=get_encoding_dtype(mz_encoding))
    unique_indices, inverse_indices = np.unique(np.digitize(mz_array, bins), return_inverse=True)
    bin_counts = np.bincount(inverse_indices)
    np.place(bin_counts, bin_counts < 1, [1])
    mz_array = np.bincount(inverse_indices, weights=mz_array) / bin_counts
    intensity_array = np.bincount(inverse_indices, weights=intensity_array)
    return mz_array, intensity_array