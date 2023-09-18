from pyTDFSDK.ctypes_data_structures import *


def tims_vis_calculate_async(tdf_sdk, handle, extraction_filter, heatmap_sizes):
    """
    Start asyncrhonous calculations for the specified data filter and heatmap sizes. Cancels any currently running
    calculations. You can get intermediate line plots and heatmaps at any time during computation. Use
    pyTDFSDK.tims_vis.tims_vis_get_state(), pyTDFSDK.tims_vis.tims_vis_wait(), or
    pyTDFSDK.tims_vis.tims_vis_wait_complete() to check/wait for the computation to be done.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param extraction_filter: The filter to be used for visualization.
    :type extraction_filter: pyTDFSDK.ctypes_data_structures.TimsVisExtractionFilter
    :param heatmap_sizes: The size definition of the 3 heatmaps.
    :type heatmap_sizes: pyTDFSDK.ctypes_data_structures.TimsVisHeatmapSizes
    :return: 1 on success, 0 on error.
    :rtype: int
    """
    try:
        result = tdf_sdk.tims_vis_calculate_async(handle, extraction_filter, heatmap_sizes)
        if result == 1:
            print('Calculation successful.')
        elif result == 0:
            print('Calculation failed.')
    except:
        raise RuntimeWarning('Unable to perform calculation.')


def tims_vis_cancel(tdf_sdk, handle):
    """
    Cancel the running calculation.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    """
    tdf_sdk.tims_vis_cancel(handle)


def tims_vis_close(tdf_sdk, handle, conn):
    """
    Close TDF/TSF dataset.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param conn: SQL database connection to analysis.tdf.
    :type conn: sqlite3.Connection
    :return: Tuple of the handle and connection.
    :rtype: tuple
    """
    if handle is not None:
        tdf_sdk.tims_vis_close(handle)
        handle = None
    if conn is not None:
        conn.close()
        conn = None
    return handle, conn


def tims_vis_get_chromatogram_line_plot(tdf_sdk, handle, x, y, line_array, transformation, size):
    """
    Get the last computed chromatogram line plot (may be an intermediate result). If the provided array is too small
    for the result data, no data is written.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param x: The requested X dimension of the plot in pixels.
    :type x: int
    :param y: The requested Y dimension of the plot in pixels.
    :type y: int
    :param line_array: The resulting array of pyTDFSDK.ctypes_data_structures.TimsVisLine.
    :type line_array: pyTDFSDK.ctypes_data_structures.TimsVisLine
    :param transformation: A transformation that might be applied to the plot.
    :type transformation: enum.Enum
    :param size: Size of the provided array.
    :type size: int
    :return: The required size of the line array or 0 on error.
    :rtype: int
    """
    pass


def tims_vis_get_image_mob_mz(tdf_sdk,
                              handle,
                              image_array,
                              transformation,
                              size):
    """
    Get the last computed (may be an intermediate result) heatmap with 1/K0 on the X axis and m/z on the Y axis. The
    dimension of the heatmap is defined in the pyTDFSDK.tims_vis.tims_vis_calculate_async() call. If the provided array
    is too small for the result data, no data is written.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param image_array: The resulting array of pyTDFSDK.ctypes_data_structures.TimsVisLine
    :type image_array: pyTDFSDK.ctypes_data_structures.TimsVisLine
    :param transformation: A transformation that might be applied to the plot.
    :type transformation: enum.Enum
    :param size: Size of the provided array.
    :type size: int
    :return: The required size of the image array or 0 on error.
    :rtype: int
    """
    pass


def tims_vis_get_image_rt_mob(tdf_sdk, handle, image_array, transformation, size):
    """
    Get the last computed (may be an intermediate result) heatmap with retention time in seconds on the X axis and 1/K0
    on the Y axis. The dimension of the heatmap is defined in the pyTDFSDK.tims_vis.tims_vis_calculate_async() call. If
    the provided array is too small for the result data, no data is written.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param image_array: The resulting array of pyTDFSDK.ctypes_data_structures.TimsVisLine
    :type image_array: pyTDFSDK.ctypes_data_structures.TimsVisLine
    :param transformation: A transformation that might be applied to the plot.
    :type transformation: enum.Enum
    :param size: Size of the provided array.
    :type size: int
    :return: The required size of the image array or 0 on error.
    :rtype: int
    """
    pass


def tims_vis_get_image_rt_mz(tdf_sdk, handle, image_array, transformation, size):
    """
    Get the last computed (may be an intermediate result) heatmap with retention time in seconds on the X axis and m/z
    on the Y axis. The dimension of the heatmap is defined in the pyTDFSDK.tims_vis.tims_vis_calculate_async() call. If
    the provided array is too small for the result data, no data is written.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param image_array: The resulting array of pyTDFSDK.ctypes_data_structures.TimsVisLine
    :type image_array: pyTDFSDK.ctypes_data_structures.TimsVisLine
    :param transformation: A transformation that might be applied to the plot.
    :type transformation: enum.Enum
    :param size: Size of the provided array.
    :type size: int
    :return: The required size of the image array or 0 on error.
    :rtype: int
    """
    pass


def tims_vis_get_mobilogram_line_plot(tdf_sdk, handle, x, y, line_array, transformation, size):
    """
    Get the last computed mobilogram line plot (may be an intermediate result). If the provided array is too small for
    the result data, no data is written.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param x: The requested X dimension of the plot in pixels.
    :type x: int
    :param y: The requested Y dimension of the plot in pixels.
    :type y: int
    :param line_array: The resulting array of pyTDFSDK.ctypes_data_structures.TimsVisLine.
    :type line_array: pyTDFSDK.ctypes_data_structures.TimsVisLine
    :param transformation: A transformation that might be applied to the plot.
    :type transformation: enum.Enum
    :param size: Size of the provided array.
    :type size: int
    :return: The required size of the line array or 0 on error.
    :rtype: int
    """
    pass


def tims_vis_get_spectrum_line_plot(tdf_sdk, handle, x, y, line_array, transformation, size):
    """
    Get the last computed spectrum line plot (may be an intermediate result). If the provided array is too small for
    the result data, no data is written.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param x: The requested X dimension of the plot in pixels.
    :type x: int
    :param y: The requested Y dimension of the plot in pixels.
    :type y: int
    :param line_array: The resulting array of pyTDFSDK.ctypes_data_structures.TimsVisLine.
    :type line_array: pyTDFSDK.ctypes_data_structures.TimsVisLine
    :param transformation: A transformation that might be applied to the plot.
    :type transformation: enum.Enum
    :param size: Size of the provided array.
    :type size: int
    :return: The required size of the line array or 0 on error.
    :rtype: int
    """
    pass


def tims_vis_get_state(tdf_sdk, handle, job_id, progress, complete):
    """
    Get the state of the last computation.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param job_id: The job ID of the last computation (an increasing counter).
    :type job_id: int
    :param progress: The progress of the last job.
    :param complete: Whether or not the last job has been completed.
    """
    pass


def tims_vis_open(tdf_sdk, bruker_d_folder_name, use_recalibrated_state=True):
    """
    Open TDF/TSF dataset and return a non-zero instance handle to be passed to subsequent API calls.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param bruker_d_folder_name: Path to a Bruker .d file containing analysis.tdf/analysis.tsf.
    :type bruker_d_folder_name: str
    :param use_recalibrated_state: Whether to use recalibrated data (True) or not (False), defaults to True.
    :type use_recalibrated_state: bool
    :return: Non-zero instance handle.
    :rtype: int
    """
    handle = tdf_sdk.tims_vis_open(bruker_d_folder_name.encode('utf-8'),
                                   1 if use_recalibrated_state else 0)
    return handle


def tims_vis_wait(tdf_sdk, handle, time_in_ms):
    """
    Wait until the current job is finished or the timeout in milliseconds is elapsed.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    :param time_in_ms: The maximum wait time in milliseconds.
    :type time_in_ms: int
    """
    tdf_sdk.tims_vis_wait(handle, time_in_ms)


def tims_vis_wait_complete(tdf_sdk, handle):
    """
    Wait until the current job is finished.

    :param tdf_sdk: Instance of TDF-SDK.
    :type tdf_sdk: ctypes.CDLL
    :param handle: Handle value for TDF/TSF dataset initialized using pyTDFSDK.tims.tims_vis_open().
    :type handle: int
    """
    tdf_sdk.tims_vis_wait_complete(handle)
