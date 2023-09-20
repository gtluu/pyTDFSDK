import numpy as np
from ctypes import POINTER, c_double, c_float, c_uint32
from pyTDFSDK.error import throw_last_tsfdata_error
from pyTDFSDK.util import call_conversion_func


def tsf_close(tdf_sdk, handle, conn):
    """
    Close TSF dataset.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param conn: SQL database connection to analysis.tsf.
    :type conn: sqlite3.Connection
    :return: Tuple of the handle and connection.
    :rtype: tuple
    """
    if handle is not None:
        tdf_sdk.tsf_close(handle)
        handle = None
    if conn is not None:
        conn.close()
        conn = None
    return handle, conn


def tsf_has_recalibrated_state(tdf_sdk, handle):
    """
    Check if the raw data has been recalibrated after acquisition (i.e. in the Bruker DataAnalysis software). Note that
    masses and 1/K0 values in the raw data SQLite files are always in the raw calibration state, not the recalibrated
    state.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :return: 1 if the raw data has been recalibrated after acquisition or 0 if not.
    :rtype: int
    """
    return tdf_sdk.tsf_has_recalibrated_state(handle)


def tsf_index_to_mz(tdf_sdk, handle, frame_id, indices):
    """
    Convert (possibly non-integer) index values for the mass dimension to m/z values.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param indices: Array of index values to be converted to m/z values.
    :type indices: numpy.array
    :return: Array of m/z values.
    :rtype: numpy.array
    """
    func = tdf_sdk.tsf_index_to_mz
    return call_conversion_func(tdf_sdk, handle, frame_id, indices, func)


def tsf_mz_to_index(tdf_sdk, handle, frame_id, mzs):
    """
    Convert m/z values to (possibly non-integer) index values for the mass dimension.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param mzs: Array of m/z values to be converted to index values.
    :type mzs: numpy.array
    :return: Array of index values.
    :rtype: numpy.array
    """
    func = tdf_sdk.tsf_mz_to_index
    return call_conversion_func(tdf_sdk, handle, frame_id, mzs, func)


def tsf_open(tdf_sdk, bruker_d_folder_name, use_recalibrated_state=True):
    """
    Open TSF dataset and return a non-zero instance handle to be passed to subsequent API calls.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param bruker_d_folder_name: Path to a Bruker .d file containing analysis.tsf.
    :type bruker_d_folder_name: str
    :param use_recalibrated_state: Whether to use recalibrated data (True) or not (False), defaults to True.
    :type use_recalibrated_state: bool
    :return: Non-zero instance handle.
    :rtype: int
    """
    handle = tdf_sdk.tsf_open(bruker_d_folder_name.encode('utf-8'), 1 if use_recalibrated_state else 0)
    return handle


def tsf_read_line_spectrum(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    """
    Read peak picked spectra for a frame from a non-TIMS (TSF) dataset.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param profile_buffer_size: Initial number of buffer bytes necessary for the output, defaults to 1024.
    :type profile_buffer_size: int
    :return: Tuple containing an array of m/z values and an array of detector counts or -1 on error.
    :rtype: tuple[numpy.array] | int
    """
    while True:
        cnt = int(profile_buffer_size)
        index_buf = np.empty(shape=cnt, dtype=np.float64)
        intensity_buf = np.empty(shape=cnt, dtype=np.float32)
        required_len = tdf_sdk.tsf_read_line_spectrum(handle,
                                                      frame_id,
                                                      index_buf.ctypes.data_as(POINTER(c_double)),
                                                      intensity_buf.ctypes.data_as(POINTER(c_float)),
                                                      profile_buffer_size)
        if required_len > profile_buffer_size:
            if required_len > 16777216:
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break
    return index_buf[0:required_len], intensity_buf[0:required_len]


def tsf_read_line_spectrum_v2(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    """
    Read peak picked spectra for a frame from a non-TIMS (TSF) dataset.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param profile_buffer_size: Initial number of buffer bytes necessary for the output, defaults to 1024.
    :type profile_buffer_size: int
    :return: Tuple containing an array of m/z values and an array of detector counts or -1 on error.
    :rtype: tuple[numpy.array] | int
    """
    while True:
        cnt = int(profile_buffer_size)
        index_buf = np.empty(shape=cnt, dtype=np.float64)
        intensity_buf = np.empty(shape=cnt, dtype=np.float32)
        required_len = tdf_sdk.tsf_read_line_spectrum_v2(handle,
                                                         frame_id,
                                                         index_buf.ctypes.data_as(POINTER(c_double)),
                                                         intensity_buf.ctypes.data_as(POINTER(c_float)),
                                                         profile_buffer_size)
        if required_len < 0:
            throw_last_tsfdata_error(tdf_sdk)
        if required_len > profile_buffer_size:
            if required_len > 16777216:
                # arbitrary limit for now...
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len  # grow buffer
        else:
            break
    return index_buf[0:required_len], intensity_buf[0:required_len]


def tsf_read_line_spectrum_with_width_v2(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    """
    Read peak picked spectra with peak width for a frame from a non-TIMS (TSF) dataset.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param profile_buffer_size: Initial number of buffer bytes necessary for the output, defaults to 1024.
    :type profile_buffer_size: int
    :return: Tuple containing an array of m/z values, an array of detector counts, and an array of widths or -1 on
        error.
    :rtype: tuple[numpy.array] | int
    """
    while True:
        cnt = int(profile_buffer_size)
        index_buf = np.empty(shape=cnt, dtype=np.float64)
        intensity_buf = np.empty(shape=cnt, dtype=np.float32)
        width_buf = np.empty(shape=cnt, dtype=np.float32)
        required_len = tdf_sdk.tsf_read_line_spectrum_with_width_v2(handle,
                                                                    frame_id,
                                                                    index_buf.ctypes.data_as(POINTER(c_double)),
                                                                    intensity_buf.ctypes.data_as(POINTER(c_float)),
                                                                    width_buf.ctypes.data_as(POINTER(c_float)),
                                                                    profile_buffer_size)
        if required_len < 0:
            throw_last_tsfdata_error(tdf_sdk)
        if required_len > profile_buffer_size:
            if required_len > 16777216:
                # arbitrary limit for now...
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break
    return index_buf[0:required_len], intensity_buf[0:required_len], width_buf[0:required_len]


def tsf_read_profile_spectrum(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    """
    Read profile spectra for a frame from a non-TIMS (TSF) dataset.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param profile_buffer_size: Initial number of buffer bytes necessary for the output, defaults to 1024.
    :type profile_buffer_size: int
    :return: Tuple containing an array of m/z values and an array of detector counts or -1 on error.
    :rtype: tuple[numpy.array] | int
    """
    while True:
        cnt = int(profile_buffer_size)
        intensity_buf = np.empty(shape=cnt, dtype=np.uint32)
        required_len = tdf_sdk.tsf_read_profile_spectrum(handle,
                                                         frame_id,
                                                         intensity_buf.ctypes.data_as(POINTER(c_uint32)),
                                                         profile_buffer_size)
        if required_len > profile_buffer_size:
            if required_len > 16777216:
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break
    index_buf = np.arange(0, intensity_buf.size, dtype=np.float64)
    return index_buf[0:required_len], intensity_buf[0:required_len]


def tsf_read_profile_spectrum_v2(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    """
    Read profile spectra for a frame from a non-TIMS (TSF) dataset.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TSF dataset initialized using pyTDFSDK.tsf.tsf_open().
    :type handle: int
    :param frame_id: ID of the frame of interest.
    :type frame_id: int
    :param profile_buffer_size: Initial number of buffer bytes necessary for the output, defaults to 1024.
    :type profile_buffer_size: int
    :return: Tuple containing an array of m/z values and an array of detector counts or -1 on error.
    :rtype: tuple[numpy.array] | int
    """
    while True:
        cnt = int(profile_buffer_size)
        intensity_buf = np.empty(shape=cnt, dtype=np.uint32)

        required_len = tdf_sdk.tsf_read_profile_spectrum_v2(handle,
                                                            frame_id,
                                                            intensity_buf.ctypes.data_as(POINTER(c_uint32)),
                                                            profile_buffer_size)
        if required_len < 0:
            throw_last_tsfdata_error(tdf_sdk)

        if required_len > profile_buffer_size:
            if required_len > 16777216:
                # arbitrary limit for now...
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break
    index_buf = np.arange(0, intensity_buf.size, dtype=np.float64)
    return index_buf[0:required_len], intensity_buf[0:required_len]


def tsf_set_num_threads(tdf_sdk, num_threads):
    """
    Set the number of threads that this DLL is allowed to use internally. The index <-> m/z transformation is
    internally parallelized using OpenMP. This call is simply forwarded to omp_set_num_threads(). This function has no
    real effect on Linux.

    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param num_threads: Number of threads to use (>= 1)
    :type num_threads: int
    """
    tdf_sdk.tsf_set_num_threads(num_threads)
