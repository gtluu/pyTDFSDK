import numpy as np
from ctypes import create_string_buffer


def throw_last_timsdata_error(tdf_sdk):
    """
    Return the last error as a string (thread-local) while throwing exception from TDF data.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :return: Runtime error.
    """
    length = tdf_sdk.tims_get_last_error_string(None, 0)
    buf = create_string_buffer(length)
    tdf_sdk.tims_get_last_error_string(buf, length)
    raise RuntimeError(buf.value)


def throw_last_timsvis_error(tdf_sdk):
    """
    Return the last error as a string (thread-local) while throwing exception from TDF/TSF data.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :return: Runtime error.
    """
    length = tdf_sdk.tims_vis_get_last_error_string(None, 0)
    buf = create_string_buffer(length)
    tdf_sdk.tims_vis_get_last_error_string(buf, length)
    raise RuntimeError(buf.value)


def throw_last_tsfdata_error(tdf_sdk):
    """
    Return the last error as a string (thread-local) while throwing exception from TSF data.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :return: Runtime error.
    """
    length = tdf_sdk.tsf_get_last_error_string(None, 0)
    buf = create_string_buffer(length)
    tdf_sdk.tsf_get_last_error_string(buf, length)
    raise RuntimeError(buf.value)
