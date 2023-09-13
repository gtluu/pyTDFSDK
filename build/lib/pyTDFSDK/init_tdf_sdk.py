import os
import platform
from pyTDFSDK.constants import *


def init_tdf_sdk_api(bruker_api_file_name=''):
    if bruker_api_file_name == '':
        if platform.system() == 'Windows':
            bruker_api_file_name = os.path.join(os.path.dirname(__file__), 'timsdata.dll')
        elif platform.system() == 'Linux':
            bruker_api_file_name = os.path.join(os.path.dirname(__file__), 'libtimsdata.so')

    tdf_sdk = ctypes.cdll.LoadLibrary(bruker_api_file_name)

    convfunc_argtypes = [ctypes.c_uint64,
                         ctypes.c_int64,
                         ctypes.POINTER(ctypes.c_double),
                         ctypes.POINTER(ctypes.c_double),
                         ctypes.c_uint32]

    # Function 1
    tdf_sdk.tims_ccs_to_oneoverk0_for_mz.argtypes = [ctypes.c_double,
                                                     ctypes.c_int32,
                                                     ctypes.c_double]
    tdf_sdk.tims_ccs_to_oneoverk0_for_mz.restype = ctypes.c_double

    # Function 2
    tdf_sdk.tims_close.argtypes = [ctypes.c_uint64]
    tdf_sdk.tims_close.restype = None

    # Function 3
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_ext.argtypes = [ctypes.c_uint64,
                                                                       ctypes.c_int64,
                                                                       ctypes.c_uint32,
                                                                       ctypes.c_uint32,
                                                                       ctypes.c_double,
                                                                       MSMS_SPECTRUM_FUNCTOR,
                                                                       ctypes.c_void_p]
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_ext.restype = ctypes.c_uint32

    # Function 4
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_v2.argtypes = [ctypes.c_uint64,
                                                                      ctypes.c_int64,
                                                                      ctypes.c_uint32,
                                                                      ctypes.c_uint32,
                                                                      MSMS_SPECTRUM_FUNCTOR,
                                                                      ctypes.c_void_p]
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_v2.restype = ctypes.c_uint32

    # Function 5
    tdf_sdk.tims_extract_chromatograms.argtypes = [ctypes.c_uint64,
                                                   CHROMATOGRAM_JOB_GENERATOR,
                                                   CHROMATOGRAM_TRACE_SINK,
                                                   ctypes.c_void_p]
    tdf_sdk.tims_extract_chromatograms.restype = ctypes.c_uint32

    # Function 6
    tdf_sdk.tims_extract_profile_for_frame.argtypes = [ctypes.c_uint64,
                                                       ctypes.c_int64,
                                                       ctypes.c_uint32,
                                                       ctypes.c_uint32,
                                                       MSMS_PROFILE_SPECTRUM_FUNCTOR,
                                                       ctypes.c_void_p]
    tdf_sdk.tims_extract_profile_for_frame.restype = ctypes.c_uint32

    # Function 7
    tdf_sdk.tims_get_last_error_string.argtypes = [ctypes.c_char_p,
                                                   ctypes.c_uint32]
    tdf_sdk.tims_get_last_error_string.restype = ctypes.c_uint32

    # Function 8
    tdf_sdk.tims_has_recalibrated_state.argtypes = [ctypes.c_uint64]
    tdf_sdk.tims_has_recalibrated_state.restype = ctypes.c_uint32

    # Function 9
    tdf_sdk.tims_index_to_mz.argtypes = convfunc_argtypes
    tdf_sdk.tims_index_to_mz.restype = ctypes.c_uint32

    # Function 10
    tdf_sdk.tims_mz_to_index.argtypes = convfunc_argtypes
    tdf_sdk.tims_mz_to_index.restype = ctypes.c_uint32

    # Function 11
    tdf_sdk.tims_oneoverk0_to_ccs_for_mz.argtypes = [ctypes.c_double,
                                                     ctypes.c_int32,
                                                     ctypes.c_double]
    tdf_sdk.tims_oneoverk0_to_ccs_for_mz.restype = ctypes.c_double

    # Function 12
    tdf_sdk.tims_oneoverk0_to_scannum.argtypes = convfunc_argtypes
    tdf_sdk.tims_oneoverk0_to_scannum.restype = ctypes.c_uint32

    # Function 13
    tdf_sdk.tims_open.argtypes = [ctypes.c_char_p,
                                  ctypes.c_uint32]
    tdf_sdk.tims_open.restype = ctypes.c_uint64

    # Function 14
    tdf_sdk.tims_open_v2.argtypes = [ctypes.c_char_p,
                                     ctypes.c_uint32,
                                     ctypes.c_uint32]
    tdf_sdk.tims_open_v2.restype = ctypes.c_uint64

    # Funtion 15
    tdf_sdk.tims_read_pasef_msms.argtypes = [ctypes.c_uint64,
                                             ctypes.POINTER(ctypes.c_int64),
                                             ctypes.c_uint32,
                                             MSMS_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_msms.restype = ctypes.c_uint32

    # Function 16
    tdf_sdk.tims_read_pasef_msms_for_frame.argtypes = [ctypes.c_uint64,
                                                       ctypes.c_int64,
                                                       MSMS_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_msms_for_frame.restype = ctypes.c_uint32

    # Function 17
    tdf_sdk.tims_read_pasef_msms_for_frame_v2.argtypes = [ctypes.c_uint64,
                                                          ctypes.c_int64,
                                                          MSMS_SPECTRUM_FUNCTION,
                                                          ctypes.POINTER(ctypes.c_void_p)]
    tdf_sdk.tims_read_pasef_msms_for_frame_v2.restype = ctypes.c_uint32

    # Function 18
    tdf_sdk.tims_read_pasef_msms_v2.argtypes = [ctypes.c_uint64,
                                                ctypes.POINTER(ctypes.c_int64),
                                                ctypes.c_uint32,
                                                MSMS_SPECTRUM_FUNCTION,
                                                ctypes.POINTER(ctypes.c_void_p)]
    tdf_sdk.tims_read_pasef_msms_v2.restype = ctypes.c_uint32

    # Function 19
    tdf_sdk.tims_read_pasef_profile_msms.argtypes = [ctypes.c_uint64,
                                                     ctypes.POINTER(ctypes.c_int64),
                                                     ctypes.c_uint32,
                                                     MSMS_PROFILE_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_profile_msms.restype = ctypes.c_uint32

    # Function 20
    tdf_sdk.tims_read_pasef_profile_msms_for_frame.argtypes = [ctypes.c_uint64,
                                                               ctypes.c_int64,
                                                               MSMS_PROFILE_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_profile_msms_for_frame.restype = ctypes.c_uint32

    # Function 21
    tdf_sdk.tims_read_pasef_profile_msms_for_frame_v2.argtypes = [ctypes.c_uint64,
                                                                  ctypes.c_int64,
                                                                  MSMS_PROFILE_SPECTRUM_FUNCTION,
                                                                  ctypes.POINTER(ctypes.c_void_p)]
    tdf_sdk.tims_read_pasef_profile_msms_for_frame_v2.restype = ctypes.c_uint32

    # Function 22
    tdf_sdk.tims_read_pasef_profile_msms_v2.argtypes = [ctypes.c_uint64,
                                                        ctypes.POINTER(ctypes.c_int64),
                                                        ctypes.c_uint32,
                                                        MSMS_PROFILE_SPECTRUM_FUNCTION,
                                                        ctypes.POINTER(ctypes.c_void_p)]
    tdf_sdk.tims_read_pasef_profile_msms_v2.restype = ctypes.c_uint32

    # Function 23 - Not Available
    #tdf_sdk.tims_read_scans_internal.argtypes = None
    #tdf_sdk.tims_read_scans_internal.restype = None

    # Function 24
    tdf_sdk.tims_read_scans_v2.argtypes = [ctypes.c_uint64,
                                           ctypes.c_int64,
                                           ctypes.c_uint32,
                                           ctypes.c_uint32,
                                           ctypes.c_void_p,
                                           ctypes.c_uint32]
    tdf_sdk.tims_read_scans_v2.restype = ctypes.c_uint32

    # Function 25
    tdf_sdk.tims_scannum_to_oneoverk0.argtypes = convfunc_argtypes
    tdf_sdk.tims_scannum_to_oneoverk0.restype = ctypes.c_uint32

    # Function 26
    tdf_sdk.tims_scannum_to_voltage.argtypes = convfunc_argtypes
    tdf_sdk.tims_scannum_to_voltage.restype = ctypes.c_uint32

    # Function 27
    tdf_sdk.tims_set_num_threads.argtypes = [ctypes.c_uint32]
    tdf_sdk.tims_set_num_threads.restype = None

    # Function 28
    tdf_sdk.tims_vis_calculate_async.argtypes = [ctypes.c_uint64,
                                                 ctypes.POINTER(TimsVisExtractionFilter),
                                                 ctypes.POINTER(TimsVisHeatmapSizes)]
    tdf_sdk.tims_vis_calculate_async.restype = ctypes.c_uint64

    # Function 29
    tdf_sdk.tims_vis_cancel.argtypes = [ctypes.c_uint64]
    tdf_sdk.tims_vis_cancel.restype = ctypes.c_uint64

    # FUnction 30
    tdf_sdk.tims_vis_close.argtypes = [ctypes.c_uint64]
    tdf_sdk.tims_vis_close.restype = None

    # Function 31
    tdf_sdk.tims_vis_get_chromatogram_line_plot.argtypes = [ctypes.c_uint64,
                                                            ctypes.c_int32,
                                                            ctypes.c_int32,
                                                            ctypes.c_uint32,
                                                            ctypes.POINTER(TimsVisLine),
                                                            ctypes.c_uint32]
    tdf_sdk.tims_vis_get_chromatogram_line_plot.restype = ctypes.c_uint32

    # Function 32
    tdf_sdk.tims_vis_get_image_mob_mz.argtypes = [ctypes.c_uint64,
                                                  ctypes.c_uint32,
                                                  ctypes.POINTER(ctypes.c_float),
                                                  ctypes.c_uint32]
    tdf_sdk.tims_vis_get_image_mob_mz.restype = ctypes.c_uint32

    # Function 33
    tdf_sdk.tims_vis_get_image_rt_mob.argtypes = [ctypes.c_uint64,
                                                  ctypes.c_uint32,
                                                  ctypes.POINTER(ctypes.c_float),
                                                  ctypes.c_uint32]
    tdf_sdk.tims_vis_get_image_rt_mob.restype = ctypes.c_uint32

    # Function 34
    tdf_sdk.tims_vis_get_image_rt_mz.argtypes = [ctypes.c_uint64,
                                                 ctypes.c_uint32,
                                                 ctypes.POINTER(ctypes.c_float),
                                                 ctypes.c_uint32]
    tdf_sdk.tims_vis_get_image_rt_mz.restype = ctypes.c_uint32

    # Function 35
    tdf_sdk.tims_vis_get_last_error_string.argtypes = [ctypes.c_char_p,
                                                       ctypes.c_uint32]
    tdf_sdk.tims_vis_get_last_error_string.restype = ctypes.c_uint32

    # Function 36
    tdf_sdk.tims_vis_get_mobilogram_line_plot.argtypes = [ctypes.c_uint64,
                                                          ctypes.c_int32,
                                                          ctypes.c_int32,
                                                          ctypes.c_uint32,
                                                          ctypes.POINTER(TimsVisLine),
                                                          ctypes.c_uint32]
    tdf_sdk.tims_vis_get_mobilogram_line_plot.restype = ctypes.c_uint32

    # Function 37
    tdf_sdk.tims_vis_get_spectrum_line_plot.argtypes = [ctypes.c_uint64,
                                                        ctypes.c_int32,
                                                        ctypes.c_int32,
                                                        ctypes.c_uint32,
                                                        ctypes.POINTER(TimsVisLine),
                                                        ctypes.c_uint32]
    tdf_sdk.tims_vis_get_spectrum_line_plot.restype = ctypes.c_uint32

    # Function 38
    tdf_sdk.tims_vis_get_state.argtypes = [ctypes.c_uint64,
                                           ctypes.c_uint64,
                                           ctypes.c_float,
                                           ctypes.c_bool]
    tdf_sdk.tims_vis_get_state_restype = ctypes.c_uint32

    # Function 39
    tdf_sdk.tims_vis_open.argtypes = [ctypes.c_char_p,
                                      ctypes.c_uint32]
    tdf_sdk.tims_vis_open.restype = ctypes.c_uint64

    # Function 40
    tdf_sdk.tims_vis_wait.argtypes = [ctypes.c_uint64,
                                      ctypes.c_uint32]
    tdf_sdk.tims_vis_wait.restype = ctypes.c_uint32

    # Function 41
    tdf_sdk.tims_vis_wait_complete.argtypes = [ctypes.c_uint64]
    tdf_sdk.tims_vis_wait_complete.restype = ctypes.c_uint32

    # Function 42
    tdf_sdk.tims_voltage_to_scannum.argtypes = convfunc_argtypes
    tdf_sdk.tims_voltage_to_scannum.restype = ctypes.c_uint32

    # Function 43
    tdf_sdk.tsf_close.argtypes = [ctypes.c_uint64]
    tdf_sdk.tsf_close.restype = None

    # Function 44
    tdf_sdk.tsf_get_last_error_string.argtypes = [ctypes.c_char_p,
                                                  ctypes.c_uint32]
    tdf_sdk.tsf_get_last_error_string.restype = ctypes.c_uint32

    # Function 45
    tdf_sdk.tsf_has_recalibrated_state.argtypes = [ctypes.c_uint64]
    tdf_sdk.tsf_has_recalibrated_state.restype = ctypes.c_uint32

    # Function 46
    tdf_sdk.tsf_index_to_mz.argtypes = convfunc_argtypes
    tdf_sdk.tsf_index_to_mz.restype = ctypes.c_uint32

    # Function 47
    tdf_sdk.tsf_mz_to_index.argtypes = convfunc_argtypes
    tdf_sdk.tsf_mz_to_index.restype = ctypes.c_uint32

    # Function 48
    tdf_sdk.tsf_open.argtypes = [ctypes.c_char_p,
                                 ctypes.c_uint32]
    tdf_sdk.tsf_open.restype = ctypes.c_uint64

    # Function 49
    tdf_sdk.tsf_read_line_spectrum.argtypes = [ctypes.c_uint64,
                                               ctypes.c_int64,
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_float),
                                               ctypes.c_uint32]
    tdf_sdk.tsf_read_line_spectrum.restype = ctypes.c_uint32

    # Function 50
    tdf_sdk.tsf_read_line_spectrum_v2.argtypes = [ctypes.c_uint64,
                                                  ctypes.c_int64,
                                                  ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_float),
                                                  ctypes.c_int32]
    tdf_sdk.tsf_read_line_spectrum_v2.restype = ctypes.c_int32

    # Function 51 - Not Available
    #tdf_sdk.tsf_read_line_spectrum_with_width.argtypes = None
    #tdf_sdk.tsf_read_line_spectrum_with_width.restype = None

    # Function 52
    tdf_sdk.tsf_read_line_spectrum_with_width_v2.argtypes = [ctypes.c_uint64,
                                                             ctypes.c_int64,
                                                             ctypes.POINTER(ctypes.c_double),
                                                             ctypes.POINTER(ctypes.c_float),
                                                             ctypes.POINTER(ctypes.c_float),
                                                             ctypes.c_int32]
    tdf_sdk.tsf_read_line_spectrum_with_width_v2.restype = ctypes.c_int32

    # Function 53
    tdf_sdk.tsf_read_profile_spectrum.argtypes = [ctypes.c_uint64,
                                                  ctypes.c_int64,
                                                  ctypes.POINTER(ctypes.c_uint32),
                                                  ctypes.c_uint32]
    tdf_sdk.tsf_read_profile_spectrum.restype = ctypes.c_uint32

    # Function 54
    tdf_sdk.tsf_read_profile_spectrum_v2.argtypes = [ctypes.c_uint64,
                                                     ctypes.c_int64,
                                                     ctypes.POINTER(ctypes.c_uint32),
                                                     ctypes.c_int32]
    tdf_sdk.tsf_read_profile_spectrum_v2.restype = ctypes.c_int32

    # Function 55
    tdf_sdk.tsf_set_num_threads.argtypes = [ctypes.c_uint32]
    tdf_sdk.tsf_set_num_threads.restype = None

    return tdf_sdk
