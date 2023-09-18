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
