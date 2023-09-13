import numpy as np
from pyTDFSDK.constants import *


# Method 7
'''
    /// Return the last error as a string (thread-local).
    ///
    /// param buf pointer to a buffer into which the error string will be written.
    ///
    /// param length length of the buffer
    ///
    /// returns the actual length of the error message (including the final zero
    /// byte). If this is longer than the input parameter 'length', you know that the
    /// returned error string was truncated to fit in the provided buffer.
'''
# Throw last TimsData error string as an exception.
def _throw_last_timsdata_error(tdf_sdk):
    length = tdf_sdk.tims_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(length)
    tdf_sdk.tims_get_last_error_string(buf, length)
    raise RuntimeError(buf.value)


# Method 35
'''
    /// Return the last error as a string (thread-local).
    ///
    /// param buf pointer to a buffer into which the error string will be written.
    ///
    /// param len length of the buffer
    ///
    /// returns the actual length of the error message (including the final zero
    /// byte). If this is longer than the input parameter 'len', you know that the
    /// returned error string was truncated to fit in the provided buffer.
'''
def _throw_last_timsvis_error(tdf_sdk):
    length = tdf_sdk.tims_vis_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(length)
    tdf_sdk.tims_vis_get_last_error_string(buf, length)
    raise RuntimeError(buf.value)


# Method 44
'''
    /// Return the last error as a string (thread-local).
    ///
    /// param buf pointer to a buffer into which the error string will be written.
    ///
    /// param len length of the buffer
    ///
    /// returns the actual length of the error message (including the final zero
    /// byte). If this is longer than the input parameter 'len', you know that the
    /// returned error string was truncated to fit in the provided buffer.
'''
def _throw_last_tsfdata_error(tdf_sdk):
    length = tdf_sdk.tsf_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(length)
    tdf_sdk.tsf_get_last_error_string(buf, length)
    raise RuntimeError(buf.value)


def __call_conversion_func(tdf_sdk, handle, frame_id, input_data, func):
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
                   in_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                   out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                   cnt)

    if success == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return out


# Method 1
'''
    /// Converts the CCS (in Angstrom^2) to 1/K0 using the Mason-Shamp equation
    ///
    /// param ccs the ccs value in Angstrom^2
    /// param charge the charge
    /// param mz the mz of the ion
    ///
    /// returns the 1/K0 value in Vs/cm2
'''
def tims_ccs_to_oneoverk0_for_mz(tdf_sdk, ccs, charge, mz):
    return tdf_sdk.tims_ccs_to_oneoverk0_for_mz(ccs, charge, mz)


# Method 2
'''
    /// Close data set.
    ///
    /// param handle obtained by tims_open(); passing 0 is ok and has no effect.
'''
def tims_close(tdf_sdk, handle, conn):
    if handle is not None:
        tdf_sdk.tims_close(handle)
        handle = None
    if conn is not None:
        conn.close()
        conn = None
    return handle, conn


# Method 3
'''
    /// Read peak-picked spectra for a tims frame with a custom peak picker resolution.
    ///
    /// Same as tims_extract_centroided_spectrum_for_frame_v2(), 
    /// but a user supplied resolution for the peak picker is applied.
    /// Can be used to prevent invalid split peaks in case of low ion statistics.
    /// The default suggested value in tims_extract_centroided_spectrum_for_frame_v2()
    /// is determined by the GlobalMetadata entry "PeakWidthEstimateValue" as
    /// 1 / PeakWidthEstimateValue for "PeakWidthEstimateType" = 1.
    ///
    /// Note: Result callback identical to the tims_read_pasef_msms_v2 methods, but
    /// only returns a single result and the parameter id is the frame_id
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_extract_centroided_spectrum_for_frame_ext(tdf_sdk, handle, frame_id, scan_begin, scan_end,
                                                   peak_picker_resolution):
    result = None

    @MSMS_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
        nonlocal result
        result = (mz_values[0:num_peaks], area_values[0:num_peaks])

    rc = tdf_sdk.tims_extract_centroided_spectrum_for_frame_ext(handle,
                                                                frame_id,
                                                                scan_begin,
                                                                scan_end,
                                                                peak_picker_resolution,
                                                                callback_for_dll,
                                                                None)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 4
'''
    /// Read peak-picked spectra for a tims frame.
    ///
    /// Given a frame ID, this function reads the frame, 
    /// sums up the corresponding scan-number ranges into a synthetic profile
    /// spectrum, performs centroiding using an algorithm and parameters
    /// suggested by Bruker, and returns the resulting spectrum (exactly one for 
    /// the frame ID).
    ///
    /// Note: Result callback identical to the tims_read_pasef_msms_v2 methods, but
    /// only returns a single result and the parameter id is the frame_id
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_extract_centroided_spectrum_for_frame_v2(tdf_sdk, handle, frame_id, scan_begin, scan_end):
    result = None

    @MSMS_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
        nonlocal result
        result = (mz_values[0:num_peaks], area_values[0:num_peaks])

    rc = tdf_sdk.tims_extract_centroided_spectrum_for_frame_v2(handle,
                                                               frame_id,
                                                               scan_begin,
                                                               scan_end,
                                                               callback_for_dll,
                                                               None)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 5
'''
    /// Extract several (MS1-only) chromatograms from an analysis.
    ///
    /// The DLL retrieves the jobs (i.e., the chromatogram definitions) from the specified generator
    /// function while iterating through the analysis. The jobs must be delivered in the order of
    /// ascending 'time_begin'.
    ///
    /// The DLL delivers chromatogram traces to the specified sink callback as soon as they are
    /// finished. When an error occurs, some of the jobs "pulled" so far might not be answered.
    ///
    /// returns 0 on error
'''
def tims_extract_chromatograms(tdf_sdk, handle, jobs, trace_sink):
    @CHROMATOGRAM_JOB_GENERATOR
    def wrap_gen(job, user_data):
        try:
            job[0] = next(jobs)
            return 1
        except StopIteration:
            return 2
        except Exception as e:
            print('extract_chromatograms: generator produced exception ', e)
            return 0

    @CHROMATOGRAM_TRACE_SINK
    def wrap_sink(job_id, num_points, frame_ids, values, user_data):
        try:
            trace_sink(job_id,
                       np.array(frame_ids[0:num_points], dtype=np.int64),
                       np.array(values[0:num_points], dtype=np.uint64))
            return 1
        except Exception as e:
            print('extract_chromatograms: sink produced exception ', e)
            return 0

    unused_user_data = 0
    rc = tdf_sdk.tims_extract_chromatogram(handle, wrap_gen, wrap_sink, unused_user_data)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)


# Method 6
'''
    /// Read "quasi profile" spectra for a tims frame.
    ///
    /// Given a frame ID, this function reads the frame, 
    /// and sums up the corresponding scan-number ranges into a synthetic profile
    /// spectrum. These "quasi" profile spectrum is passed back.
    ///
    /// Note: Result callback identical to the tims_read_pasef_profile_msms_v2 methods,
    /// but only returns a single result and the parameter id is the frame_id
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_extract_profile_for_frame(tdf_sdk, handle, frame_id, scan_begin, scan_end):
    result = None

    @MSMS_PROFILE_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_points, intensity_values):
        nonlocal result
        result = intensity_values[0:num_points]

    rc = tdf_sdk.tims_extract_profile_for_frame(handle, frame_id, scan_begin, scan_end, callback_for_dll, None)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 8
'''
    /// Returns 1 if the raw data have been recalibrated after acquisition, e.g. in the
    /// DataAnalysis software. Note that masses and 1/K0 values in the raw-data SQLite
    /// file are always in the raw calibration state, not the recalibrated state.
'''
def tims_has_recalibrated_state(tdf_sdk, handle):
    return tdf_sdk.tims_has_recalibrated_state(handle)


# Method 9
'''
    /// m/z transformation: convert back and forth between (possibly non-integer) index
    /// values and m/z values.
'''
def tims_index_to_mz(tdf_sdk, handle, frame_id, indices):
    func = tdf_sdk.tims_index_to_mz
    return __call_conversion_func(tdf_sdk, handle, frame_id, indices, func)


# Method 10
'''
    /// m/z transformation: convert back and forth between (possibly non-integer) index
    /// values and m/z values.
'''
def tims_mz_to_index(tdf_sdk, handle, frame_id, mzs):
    func = tdf_sdk.tims_mz_to_index
    return __call_conversion_func(tdf_sdk, handle, frame_id, mzs, func)


# Method 11
'''
    /// Converts the 1/K0 value to CCS (in Angstrom^2) using the Mason-Shamp equation

    /// param ook0 the 1/K0 value in Vs/cm2
    /// param charge the charge
    /// param mz the mz of the ion

    /// returns the CCS value in Angstrom^2
'''
def tims_oneoverk0_to_ccs_for_mz(tdf_sdk, ook0, charge, mz):
    return tdf_sdk.tims_oneoverk0_to_ccs_for_mz(ook0, charge, mz)


# Method 12
'''
    /// mobility transformation: convert back and forth between (possibly non-integer)
    /// scan numbers and 1/K0 values.
'''
def tims_oneoverk0_to_scannum(tdf_sdk, handle, frame_id, mobilities):
    func = tdf_sdk.tims_oneoverk0_to_scannum
    return __call_conversion_func(tdf_sdk, handle, frame_id, mobilities, func)


# Method 13
'''
    /// Open data set.
    ///
    /// On success, returns a non-zero instance handle that needs to be passed to
    /// subsequent API calls, in particular to the required call to tims_close().
    ///
    /// On failure, returns 0, and you can use tims_get_last_error_string() to obtain a
    /// string describing the problem.
    ///
    /// Uses NoPressureCompensation.
    ///
    /// param analysis_directory_name the name of the directory in the file system that
    /// contains the analysis data, in UTF-8 encoding.
    ///
    /// param use_recalibrated_state if non-zero, use the most recent recalibrated state
    /// of the analysis, if there is one; if zero, use the original "raw" calibration
    /// written during acquisition time.
'''
def tims_open(tdf_sdk, bruker_d_folder_name, use_recalibrated_state=True):
    handle = tdf_sdk.tims_open(bruker_d_folder_name.encode('utf-8'), 1 if use_recalibrated_state else 0)
    return handle


# Method 14
'''
    /// Open data set.
    ///
    /// On success, returns a non-zero instance handle that needs to be passed to
    /// subsequent API calls, in particular to the required call to tims_close().
    ///
    /// On failure, returns 0, and you can use tims_get_last_error_string() to obtain a
    /// string describing the problem.
    ///
    /// param analysis_directory_name the name of the directory in the file system that
    /// contains the analysis data, in UTF-8 encoding.
    ///
    /// param use_recalibrated_state if non-zero, use the most recent recalibrated state
    /// of the analysis, if there is one; if zero, use the original "raw" calibration
    /// written during acquisition time.
    ///
    /// param pressure_compensation_strategy the pressure compensation strategy
'''
def tims_open_v2(tdf_sdk, bruker_d_folder_name, pressure_compensation_strategy, use_recalibrated_state=True):
    handle = tdf_sdk.tims_open_v2(bruker_d_folder_name.encode('utf-8'),
                                  1 if use_recalibrated_state else 0,
                                  pressure_compensation_strategy.value)
    return handle


# Method 15
'''
    /// Read peak-picked MS/MS spectra for a list of PASEF precursors.
    ///
    /// Given a list of PASEF precursor IDs, this function reads all necessary PASEF
    /// frames, sums up the corresponding scan-number ranges into synthetic profile
    /// spectra for each precursor, performs centroiding using an algorithm and parameters
    /// suggested by Bruker, and returns the resulting MS/MS spectra (one for each
    /// precursor ID).
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_msms(tdf_sdk, handle, precursor_list):
    precursors_for_dll = np.array(precursor_list, dtype=np.int64)

    result = {}

    @MSMS_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
        result[precursor_id] = (mz_values[0:num_peaks], area_values[0:num_peaks])

    rc = tdf_sdk.tims_read_pasef_msms(handle,
                                      precursors_for_dll.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
                                      len(precursor_list),
                                      callback_for_dll)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 16
'''
    /// Read peak-picked MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_msms.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_msms_for_frame(tdf_sdk, handle, frame_id):
    result = {}

    @MSMS_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
        result[precursor_id] = (mz_values[0:num_peaks], area_values[0:num_peaks])

    rc = tdf_sdk.tims_read_pasef_msms_for_frame(handle, frame_id, callback_for_dll)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 17
'''
    /// Read peak-picked MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_msms.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_msms_for_frame_v2(tdf_sdk, handle, frame_id):
    result = {}

    @MSMS_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
        nonlocal result
        result[precursor_id] = (mz_values[0:num_peaks], area_values[0:num_peaks])

    rc = tdf_sdk.tims_read_pasef_msms_for_frame_v2(handle, frame_id, callback_for_dll, None)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 18
'''
    /// Read peak-picked MS/MS spectra for a list of PASEF precursors.
    ///
    /// Given a list of PASEF precursor IDs, this function reads all necessary PASEF
    /// frames, sums up the corresponding scan-number ranges into synthetic profile
    /// spectra for each precursor, performs centroiding using an algorithm and parameters
    /// suggested by Bruker, and returns the resulting MS/MS spectra (one for each
    /// precursor ID).
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_msms_v2(tdf_sdk, handle, precursor_list):
    precursors_for_dll = np.array(precursor_list, dtype=np.int64)

    result = {}

    @MSMS_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
        nonlocal result
        result[precursor_id] = (mz_values[0:num_peaks], area_values[0:num_peaks])

    rc = tdf_sdk.tims_read_pasef_msms_v2(handle,
                                         precursors_for_dll.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
                                         len(precursor_list),
                                         callback_for_dll,
                                         None)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 19
'''
    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a list of PASEF precursor IDs, this function reads all necessary PASEF
    /// frames, sums up the corresponding scan-number ranges into synthetic profile
    /// spectra for each precursor. These "quasi" profile spectra are passed back - one
    /// for each precursor ID.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_profile_msms(tdf_sdk, handle, precursor_list):
    precursors_for_dll = np.array(precursor_list, dtype=np.int64)

    result = {}

    @MSMS_PROFILE_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_points, intensity_values):
        result[precursor_id] = intensity_values[0:num_points]

    rc = tdf_sdk.tims_read_pasef_profile_msms(handle,
                                              precursors_for_dll.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
                                              len(precursor_list),
                                              callback_for_dll)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 20
'''
    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads for all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_profile_msms.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_profile_msms_for_frame(tdf_sdk, handle, frame_id):
    result = {}

    @MSMS_PROFILE_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_points, intensity_values):
        result[precursor_id] = intensity_values[0:num_points]

    rc = tdf_sdk.tims_read_pasef_profile_msms_for_frame(handle, frame_id, callback_for_dll)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 21
'''
    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads for all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_profile_msms.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_profile_msms_for_frame_v2(tdf_sdk, handle, frame_id):
    result = {}

    @MSMS_PROFILE_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_points, intensity_values):
        nonlocal result
        result[precursor_id] = intensity_values[0:num_points]

    rc = tdf_sdk.tims_read_pasef_profile_msms_for_frame_v2(handle, frame_id, callback_for_dll, None)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 22
'''
    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a list of PASEF precursor IDs, this function reads all necessary PASEF
    /// frames, sums up the corresponding scan-number ranges into synthetic profile
    /// spectra for each precursor. These "quasi" profile spectra are passed back - one
    /// for each precursor ID.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error
'''
def tims_read_pasef_profile_msms_v2(tdf_sdk, handle, precursor_list):
    precursors_for_dll = np.array(precursor_list, dtype=np.int64)

    result = {}

    @MSMS_PROFILE_SPECTRUM_FUNCTOR
    def callback_for_dll(precursor_id, num_points, intensity_values):
        nonlocal result
        result[precursor_id] = intensity_values[0:num_points]

    rc = tdf_sdk.tims_read_pasef_profile_msms_v2(handle,
                                                 precursors_for_dll.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
                                                 len(precursor_list),
                                                 callback_for_dll,
                                                 None)

    if rc == 0:
        _throw_last_timsdata_error(tdf_sdk)

    return result


# Method 23 - Not Available


# Method 24
'''
    /// Read a range of scans from a single frame.
    ///
    /// Output layout: (N = scan_end - scan_begin = number of requested scans)
    ///   N x uint32_t: number of peaks in each of the N requested scans
    ///   N x (two uint32_t arrays: first indices, then intensities)
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// returns 0 on error, otherwise the number of buffer bytes necessary for the output
    /// of this call (if this is larger than the provided buffer length, the result is not
    /// complete).
'''
def tims_read_scans_v2(tdf_sdk, handle, frame_id, scan_begin, scan_end, initial_frame_buffer_size=128):
    # buffer-growing loop
    while True:
        cnt = int(initial_frame_buffer_size)
        buf = np.empty(shape=cnt, dtype=np.uint32)
        length = 4 * cnt

        required_len = tdf_sdk.tims_read_scans_v2(handle,
                                                  frame_id,
                                                  scan_begin,
                                                  scan_end,
                                                  buf.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                                                  length)

        if required_len == 0:
            _throw_last_timsdata_error(tdf_sdk)

        if required_len > length:
            if required_len > 16777216:
                # arbitrary limit for now...
                raise RuntimeError('Maximum expected frame size exceeded.')
            initial_frame_buffer_size = required_len / 4 + 1  # grow buffer
        else:
            break

    result = []
    d = scan_end - scan_begin
    for i in range(scan_begin, scan_end):
        npeaks = buf[i-scan_begin]
        indices = buf[d : d+npeaks]
        d += npeaks
        intensities = buf[d : d+npeaks]
        d += npeaks
        result.append((indices,intensities))

    return result


# Method 25
'''
    /// mobility transformation: convert back and forth between (possibly non-integer)
    /// scan numbers and 1/K0 values.
'''
def tims_scannum_to_oneoverk0(tdf_sdk, handle, frame_id, scan_nums):
    func = tdf_sdk.tims_scannum_to_oneoverk0
    return __call_conversion_func(tdf_sdk, handle, frame_id, scan_nums, func)


# Method 26
'''
    /// mobility transformation: convert back and forth between (possibly non-integer)
    /// scan numbers and TIMS voltages.
'''
def tims_scannum_to_voltage(tdf_sdk, handle, frame_id, scan_nums):
    func = tdf_sdk.tims_scannum_to_voltage
    return __call_conversion_func(tdf_sdk, handle, frame_id, scan_nums, func)


# Method 27
'''
    /// Set the number of threads that this DLL is allowed to use internally. [The
    /// index<->m/z transformation is internally parallelized using OpenMP; this call is
    /// simply forwarded to omp_set_num_threads(). Has no effect on Linux].
    ///
    /// param n number of threads to use (n must be >= 1).
'''
def tims_set_num_threads(tdf_sdk, num_threads):
    tdf_sdk.tims_set_num_threads(num_threads)


# Method 28
'''
    /// Start asynchronous calculations for the specified data filter and heatmap sizes.
    /// Cancels any currently running calculations.
    /// You can get intermediate line plots and heatmaps at any time during computation.
    /// Use tims_vis_get_state, tims_vis_wait or tims_vis_wait_complete to check/wait for
    /// the computation to be done.
    ///
    /// param handle obtained by tims_open()
    ///
    /// param filter the filter use for visualization
    ///
    /// param heatmap_sizes the size definition of the 3 heatmaps
    ///
    /// returns 1 on success, 0 on error
'''
def tims_vis_calculate_async(tdf_sdk, handle, extraction_filter, heatmap_sizes):
    try:
        result = tdf_sdk.tims_vis_calculate_async(handle, extraction_filter, heatmap_sizes)
        if result == 1:
            print('Calculation successful.')
        elif result == 0:
            print('Calculation failed.')
    except:
        raise RuntimeWarning('Unable to perform calculation.')


# Method 29
'''
    /// Cancel a running calculation
    ///
    /// param handle obtained by tims_open()
'''
def tims_vis_cancel(tdf_sdk, handle):
    tdf_sdk.tims_vis_cancel(handle)


# Method 30
'''
    /// Close data set.
    ///
    /// \param handle obtained by tims_vis_open(); passing 0 is ok and has no effect.
'''
def tims_vis_close(tdf_sdk, handle, conn):
    if handle is not None:
        tdf_sdk.tims_vis_close(handle)
        handle = None
    if conn is not None:
        conn.close()
        conn = None
    return handle, conn


# Method 31
'''
    /// Get the last computed chromatogram line plot (may be an intermediate result).
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// param handle obtained by tims_open()
    /// param pixelX the requested X dimension of the plot in pixels
    /// param pixelY the requested Y dimension of the plot in pixels
    /// param transformation a transformation that might be applied to the plot
    /// param line_array the resulting array of tims_vis_line, allocated by user
    /// param size size of the provided array
    ///
    /// returns 0 on error, otherwise the required size of the line_array
'''
def tims_vis_get_chromatogram_line_plot():
    pass


# Method 32
'''
    /// Get the last computed (maybe intermediate) heatmap
    /// with 1/K0 on the X axis and mz on the Y axis.
    /// The dimension of the heatmap is defined in the tims_vis_calculate_async call.
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// param handle obtained by tims_open()
    /// param transformation a transformation that might be applied to the plot
    /// param image_array the resulting array of tims_vis_line, allocated by user
    /// param size size of the provided array
    ///
    /// returns 0 on error, otherwise the required size of the image_array
'''
def tims_vis_get_image_mob_mz(tdf_sdk, handle, image_array, transformation=TimsVisTransformation.NONE.value,
                              buffer_size=1024):
    pass


# Method 33
'''
    /// Get the last computed (maybe intermediate) heatmap
    /// with retention time in seconds on the X axis and 1/Ko on the Y axis.
    /// The dimension of the heatmap is defined in the tims_vis_calculate_async call.
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// param handle obtained by tims_open()
    /// param transformation a transformation that might be applied to the plot
    /// param image_array the resulting array of tims_vis_line, allocated by user
    /// param size size of the provided array
    ///
    /// returns 0 on error, otherwise the required size of the image_array
'''
def tims_vis_get_image_rt_mob():
    pass


# Method 34
'''
    /// Get the last computed (maybe intermediate) heatmap
    /// with retention time in seconds on the X axis and mz on the Y axis.
    /// The dimension of the heatmap is defined in the tims_vis_calculate_async call.
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// param handle obtained by tims_open()
    /// param transformation a transformation that might be applied to the plot
    /// param image_array the resulting array of tims_vis_line, allocated by user
    /// param size size of the provided array
    ///
    /// returns 0 on error, otherwise the required size of the image_array
'''
def tims_vis_get_image_rt_mz():
    pass


# Method 36
'''
    /// Get the last computed mobilogram line plot (may be an intermediate result).
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// param handle obtained by tims_open()
    /// param pixelX the requested X dimension of the plot in pixels
    /// param pixelY the requested Y dimension of the plot in pixels
    /// param transformation a transformation that might be applied to the plot
    /// param line_array the resulting array of tims_vis_line, allocated by user
    /// param size size of the provided array
    ///
    /// returns 0 on error, otherwise the required size of the line_array
'''
def tims_vis_get_mobilogram_line_plot(tdf_sdk, handle, x, y, transformation, line_array):
    pass


# Method 37
'''
    /// Get the last computed spectrum line plot (may be an intermediate result).
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// param handle obtained by tims_open()
    /// param pixelX the requested X dimension of the plot in pixels
    /// param pixelY the requested Y dimension of the plot in pixels
    /// param transformation a transformation that might be applied to the plot
    /// param line_array the resulting array of tims_vis_line, allocated by user
    /// param size size of the provided array
    ///
    /// returns 0 on error, otherwise the required size of the line_array
'''
def tims_vis_get_spectrum_line_plot():
    pass


# Method 38
'''
    /// Get the state of the last started computation.
    ///
    /// param handle obtained by tims_open()
    /// param[out] job_id the job id of the last computation (just an increasing counter)
    /// param[out] progress the progress of the last job
    /// param[out] complete if the last job has been completed
'''
def tims_vis_get_state(tdf_sdk, handle, job_id, progress, complete):
    #return tdf_sdk.tims_vis_get_state(handle, job_id, progress, complete)
    pass


# Method 39
'''
    /// Open data set.
    ///
    /// On success, returns a non-zero instance handle that needs to be passed to
    /// subsequent API calls, in particular to the required call to tims_close().
    ///
    /// On failure, returns 0, and you can use tims_vis_get_last_error_string() to obtain a
    /// string describing the problem.
    ///
    /// param analysis_directory_name the name of the directory in the file system that
    /// contains the analysis data, in UTF-8 encoding.
    ///
    /// param use_recalibrated_state if non-zero, use the most recent recalibrated state
    /// of the analysis, if there is one; if zero, use the original "raw" calibration
    /// written during acquisition time.
'''
def tims_vis_open(tdf_sdk, bruker_d_folder_name, use_recalibrated_state=True):
    handle = tdf_sdk.tims_vis_open(bruker_d_folder_name.encode('utf-8'),
                                   1 if use_recalibrated_state else 0)
    return handle


# Method 40
'''
    /// Wait until the current job is finished or the timeout in milliseconds is elapsed
    ///
    /// param handle obtained by tims_open()
    /// param tims_in_ms the maximum wait time in milliseconds
'''
def tims_vis_wait(tdf_sdk, handle, time_in_ms):
    tdf_sdk.tims_vis_wait(handle, time_in_ms)


# Method 41
'''
    /// Wait until the current job is finished
    ///
    /// param handle obtained by tims_open()
'''
def tims_vis_wait_complete(tdf_sdk, handle):
    tdf_sdk.tims_vis_wait_complete(handle)


# Method 42
'''
    /// mobility transformation: convert back and forth between (possibly non-integer)
    /// scan numbers and TIMS voltages.
'''
def tims_voltage_to_scannum(tdf_sdk, handle, frame_id, voltages):
    func = tdf_sdk.tims_voltage_to_scannum
    return __call_conversion_func(tdf_sdk, handle, frame_id, voltages, func)


# Method 43
'''
    /// Close data set.
    ///
    /// param handle obtained by tsf_open(); passing 0 is ok and has no effect.
'''
def tsf_close(tdf_sdk, handle, conn):
    if handle is not None:
        tdf_sdk.tsf_close(handle)
        handle = None
    if conn is not None:
        conn.close()
        conn = None
    return handle, conn


# Method 45
'''
    /// Returns 1 if the raw data have been recalibrated after acquisition, e.g. in the
    /// DataAnalysis software. Note that masses and 1/K0 values in the raw-data SQLite
    /// file are always in the raw calibration state, not the recalibrated state.
'''
def tsf_has_recalibrated_state(tdf_sdk, handle):
    return tdf_sdk.tsf_has_recalibrated_state(handle)


# Method 46
'''
    /// m/z transformation: convert back and forth between (possibly non-integer) index
    /// values and m/z values.
'''
def tsf_index_to_mz(tdf_sdk, handle, frame_id, indices):
    func = tdf_sdk.tsf_index_to_mz
    return __call_conversion_func(tdf_sdk, handle, frame_id, indices, func)


# Method 47
'''
    /// m/z transformation: convert back and forth between (possibly non-integer) index
    /// values and m/z values.
'''
def tsf_mz_to_index(tdf_sdk, handle, frame_id, mzs):
    func = tdf_sdk.tsf_mz_to_index
    return __call_conversion_func(tdf_sdk, handle, frame_id, mzs, func)


# Method 48
'''
    /// Open data set.
    ///
    /// On success, returns a non-zero instance handle that needs to be passed to
    /// subsequent API calls, in particular to the required call to tims_close().
    ///
    /// On failure, returns 0, and you can use tims_get_last_error_string() to obtain a
    /// string describing the problem.
    ///
    /// \param analysis_directory_name the name of the directory in the file system that
    /// contains the analysis data, in UTF-8 encoding.
    ///
    /// \param use_recalibrated_state if non-zero, use the most recent recalibrated state
    /// of the analysis, if there is one; if zero, use the original "raw" calibration
    /// written during acquisition time.
'''
def tsf_open(tdf_sdk, bruker_d_folder_name, use_recalibrated_state=True):
    handle = tdf_sdk.tims_open(bruker_d_folder_name.encode('utf-8'), 1 if use_recalibrated_state else 0)
    return handle


# Method 49
'''
    /// Read a line spectrum. Fails if no line spectrum is contained (GlobalMetatdata.HasLineSpectra == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// param[in] handle the handle used for reading
    /// param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// param[out] index_array index values as double array 
    /// param[out] intensity_array intensity values as float array 
    /// param[in] length the length of the provided arrays
    ///
    /// returns -1 on error, otherwise the number of entries necessary for the output arrays
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
'''
def tsf_read_line_spectrum(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    while True:
        cnt = int(profile_buffer_size)
        index_buf = np.empty(shape=cnt, dtype=np.float64)
        intensity_buf = np.empty(shape=cnt, dtype=np.float32)

        required_len = tdf_sdk.tsf_read_line_spectrum(handle,
                                                      frame_id,
                                                      index_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                                      intensity_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                                      profile_buffer_size)

        if required_len > profile_buffer_size:
            if required_len > 16777216:
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break

    return index_buf[0:required_len], intensity_buf[0:required_len]


# Method 50
'''
    /// Read a line spectrum. Fails if no line spectrum is contained (GlobalMetatdata.HasLineSpectra == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// param[in] handle the handle used for reading
    /// param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// param[out] index_array index values as double array 
    /// param[out] intensity_array intensity values as float array 
    /// param[in] length the length of the provided arrays
    ///
    /// returns -1 on error, otherwise the number of entries necessary for the output arrays
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
'''
def tsf_read_line_spectrum_v2(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    while True:
        cnt = int(profile_buffer_size)
        index_buf = np.empty(shape=cnt, dtype=np.float64)
        intensity_buf = np.empty(shape=cnt, dtype=np.float32)

        required_len = tdf_sdk.tsf_read_line_spectrum_v2(handle,
                                                         frame_id,
                                                         index_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                                         intensity_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                                         profile_buffer_size)

        if required_len < 0:
            _throw_last_tsfdata_error(tdf_sdk)

        if required_len > profile_buffer_size:
            if required_len > 16777216:
                # arbitrary limit for now...
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len  # grow buffer
        else:
            break

    return index_buf[0:required_len], intensity_buf[0:required_len]


# Method 51 - Not Available


# Method 52
'''
    /// Read a line spectrum. Fails if no line spectrum or no peak width is contained 
    /// (GlobalMetatdata.HasLineSpectra == 0 or GlobalMetatdata.HasLineSpectraPeakWidth == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// param[in] handle the handle used for reading
    /// param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// param[out] index_array index values as double array 
    /// param[out] intensity_array intensity values as float array 
    /// param[out] width_array width values as float array 
    /// param[in] length the length of the provided arrays
    ///
    /// returns -1 on error, otherwise the number of entries necessary for the output arrays
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
'''
def tsf_read_line_spectrum_with_width_v2(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    while True:
        cnt = int(profile_buffer_size)
        index_buf = np.empty(shape=cnt, dtype=np.float64)
        intensity_buf = np.empty(shape=cnt, dtype=np.float32)
        width_buf = np.empty(shape=cnt, dtype=np.float32)

        required_len = tdf_sdk.tsf_read_line_spectrum_with_width_v2(handle,
                                                                    frame_id,
                                                                    index_buf.ctypes.data_as(ctypes.POINTER(ctypes.double)),
                                                                    intensity_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                                                    width_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                                                    profile_buffer_size)

        if required_len < 0:
            _throw_last_tsfdata_error(tdf_sdk)

        if required_len > profile_buffer_size:
            if required_len > 16777216:
                # arbitrary limit for now...
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break

    return index_buf[0:required_len], intensity_buf[0:required_len], width_buf[0:required_len]


# Method 53
'''
    /// Read a profile spectrum. Fails if no profile is contained (GlobalMetatdata.HasProfileSpectra == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// param[in] handle the handle used for reading
    /// param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// param[out] profile_array intensity values as uint32_t array, position in the array is the index
    /// param[in] length the length of the provided array
    ///
    /// returns -1 on error, otherwise the number of entries necessary for the output array
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
'''
def tsf_read_profile_spectrum(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    while True:
        cnt = int(profile_buffer_size)
        intensity_buf = np.empty(shape=cnt, dtype=np.uint32)

        required_len = tdf_sdk.tsf_read_profile_spectrum(handle,
                                                         frame_id,
                                                         intensity_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                                                         profile_buffer_size)

        if required_len > profile_buffer_size:
            if required_len > 16777216:
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break

    index_buf = np.arange(0, intensity_buf.size, dtype=np.float64)

    return index_buf[0:required_len], intensity_buf[0:required_len]


# Method 54
'''
    /// Read a profile spectrum. Fails if no profile is contained (GlobalMetatdata.HasProfileSpectra == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// param[in] handle the handle used for reading
    /// param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// param[out] profile_array intensity values as uint32_t array, position in the array is the index
    /// param[in] length the length of the provided array
    ///
    /// returns -1 on error, otherwise the number of entries necessary for the output array
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
'''
def tsf_read_profile_spectrum_v2(tdf_sdk, handle, frame_id, profile_buffer_size=1024):
    while True:
        cnt = int(profile_buffer_size)
        intensity_buf = np.empty(shape=cnt, dtype=np.uint32)

        required_len = tdf_sdk.tsf_read_profile_spectrum_v2(handle,
                                                            frame_id,
                                                            intensity_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                                                            profile_buffer_size)

        if required_len < 0:
            _throw_last_tsfdata_error(tdf_sdk)

        if required_len > profile_buffer_size:
            if required_len > 16777216:
                # arbitrary limit for now...
                raise RuntimeError('Maximum expected frame size exceeded.')
            profile_buffer_size = required_len
        else:
            break

    return intensity_buf[0:required_len]


# Method 55
'''
    /// Set the number of threads that this DLL is allowed to use internally. [The
    /// index<->m/z transformation is internally parallelized using OpenMP; this call is
    /// simply forwarded to omp_set_num_threads()].
    ///
    /// param n number of threads to use (n must be >= 1).
'''
def tsf_set_num_threads(tdf_sdk, num_threads):
    tdf_sdk.tsf_set_num_threads(num_threads)
