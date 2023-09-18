import os
import platform
from ctypes import CFUNCTYPE, POINTER, Structure, c_int64, c_uint32, c_double, c_float, c_int32, c_void_p, c_uint64, \
    cdll, c_char_p, c_bool
from pyTDFSDK.ctypes_data_structures import *


def init_tdf_sdk_api(bruker_api_file_name=''):
    if bruker_api_file_name == '':
        if platform.system() == 'Windows':
            bruker_api_file_name = os.path.join(os.path.split(os.path.dirname(__file__))[0],
                                                'TDF-SDK',
                                                'timsdata.dll')
        elif platform.system() == 'Linux':
            bruker_api_file_name = os.path.join(os.path.split(os.path.dirname(__file__))[0],
                                                'TDF-SDK',
                                                'libtimsdata.so')

    tdf_sdk = cdll.LoadLibrary(bruker_api_file_name)

    convfunc_argtypes = [c_uint64,
                         c_int64,
                         POINTER(c_double),
                         POINTER(c_double),
                         c_uint32]

    tdf_sdk.tims_ccs_to_oneoverk0_for_mz.argtypes = [c_double,
                                                     c_int32,
                                                     c_double]
    tdf_sdk.tims_ccs_to_oneoverk0_for_mz.restype = c_double

    tdf_sdk.tims_close.argtypes = [c_uint64]
    tdf_sdk.tims_close.restype = None

    tdf_sdk.tims_extract_centroided_spectrum_for_frame_ext.argtypes = [c_uint64,
                                                                       c_int64,
                                                                       c_uint32,
                                                                       c_uint32,
                                                                       c_double,
                                                                       MSMS_SPECTRUM_FUNCTOR,
                                                                       c_void_p]
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_ext.restype = c_uint32

    tdf_sdk.tims_extract_centroided_spectrum_for_frame_v2.argtypes = [c_uint64,
                                                                      c_int64,
                                                                      c_uint32,
                                                                      c_uint32,
                                                                      MSMS_SPECTRUM_FUNCTOR,
                                                                      c_void_p]
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_v2.restype = c_uint32

    tdf_sdk.tims_extract_chromatograms.argtypes = [c_uint64,
                                                   CHROMATOGRAM_JOB_GENERATOR,
                                                   CHROMATOGRAM_TRACE_SINK,
                                                   c_void_p]
    tdf_sdk.tims_extract_chromatograms.restype = c_uint32

    tdf_sdk.tims_extract_profile_for_frame.argtypes = [c_uint64,
                                                       c_int64,
                                                       c_uint32,
                                                       c_uint32,
                                                       MSMS_PROFILE_SPECTRUM_FUNCTOR,
                                                       c_void_p]
    tdf_sdk.tims_extract_profile_for_frame.restype = c_uint32

    tdf_sdk.tims_get_last_error_string.argtypes = [c_char_p,
                                                   c_uint32]
    tdf_sdk.tims_get_last_error_string.restype = c_uint32

    tdf_sdk.tims_has_recalibrated_state.argtypes = [c_uint64]
    tdf_sdk.tims_has_recalibrated_state.restype = c_uint32

    tdf_sdk.tims_index_to_mz.argtypes = convfunc_argtypes
    tdf_sdk.tims_index_to_mz.restype = c_uint32

    tdf_sdk.tims_mz_to_index.argtypes = convfunc_argtypes
    tdf_sdk.tims_mz_to_index.restype = c_uint32

    tdf_sdk.tims_oneoverk0_to_ccs_for_mz.argtypes = [c_double,
                                                     c_int32,
                                                     c_double]
    tdf_sdk.tims_oneoverk0_to_ccs_for_mz.restype = c_double

    tdf_sdk.tims_oneoverk0_to_scannum.argtypes = convfunc_argtypes
    tdf_sdk.tims_oneoverk0_to_scannum.restype = c_uint32

    tdf_sdk.tims_open.argtypes = [c_char_p,
                                  c_uint32]
    tdf_sdk.tims_open.restype = c_uint64

    tdf_sdk.tims_open_v2.argtypes = [c_char_p,
                                     c_uint32,
                                     c_uint32]
    tdf_sdk.tims_open_v2.restype = c_uint64

    tdf_sdk.tims_read_pasef_msms.argtypes = [c_uint64,
                                             POINTER(c_int64),
                                             c_uint32,
                                             MSMS_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_msms.restype = c_uint32

    tdf_sdk.tims_read_pasef_msms_for_frame.argtypes = [c_uint64,
                                                       c_int64,
                                                       MSMS_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_msms_for_frame.restype = c_uint32

    tdf_sdk.tims_read_pasef_msms_for_frame_v2.argtypes = [c_uint64,
                                                          c_int64,
                                                          MSMS_SPECTRUM_FUNCTION,
                                                          POINTER(c_void_p)]
    tdf_sdk.tims_read_pasef_msms_for_frame_v2.restype = c_uint32

    tdf_sdk.tims_read_pasef_msms_v2.argtypes = [c_uint64,
                                                POINTER(c_int64),
                                                c_uint32,
                                                MSMS_SPECTRUM_FUNCTION,
                                                POINTER(c_void_p)]
    tdf_sdk.tims_read_pasef_msms_v2.restype = c_uint32

    tdf_sdk.tims_read_pasef_profile_msms.argtypes = [c_uint64,
                                                     POINTER(c_int64),
                                                     c_uint32,
                                                     MSMS_PROFILE_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_profile_msms.restype = c_uint32

    tdf_sdk.tims_read_pasef_profile_msms_for_frame.argtypes = [c_uint64,
                                                               c_int64,
                                                               MSMS_PROFILE_SPECTRUM_FUNCTOR]
    tdf_sdk.tims_read_pasef_profile_msms_for_frame.restype = c_uint32

    tdf_sdk.tims_read_pasef_profile_msms_for_frame_v2.argtypes = [c_uint64,
                                                                  c_int64,
                                                                  MSMS_PROFILE_SPECTRUM_FUNCTION,
                                                                  POINTER(c_void_p)]
    tdf_sdk.tims_read_pasef_profile_msms_for_frame_v2.restype = c_uint32

    tdf_sdk.tims_read_pasef_profile_msms_v2.argtypes = [c_uint64,
                                                        POINTER(c_int64),
                                                        c_uint32,
                                                        MSMS_PROFILE_SPECTRUM_FUNCTION,
                                                        POINTER(c_void_p)]
    tdf_sdk.tims_read_pasef_profile_msms_v2.restype = c_uint32

    tdf_sdk.tims_read_scans_v2.argtypes = [c_uint64,
                                           c_int64,
                                           c_uint32,
                                           c_uint32,
                                           c_void_p,
                                           c_uint32]
    tdf_sdk.tims_read_scans_v2.restype = c_uint32

    tdf_sdk.tims_scannum_to_oneoverk0.argtypes = convfunc_argtypes
    tdf_sdk.tims_scannum_to_oneoverk0.restype = c_uint32

    tdf_sdk.tims_scannum_to_voltage.argtypes = convfunc_argtypes
    tdf_sdk.tims_scannum_to_voltage.restype = c_uint32

    # Function 27
    tdf_sdk.tims_set_num_threads.argtypes = [c_uint32]
    tdf_sdk.tims_set_num_threads.restype = None

    tdf_sdk.tims_vis_calculate_async.argtypes = [c_uint64,
                                                 POINTER(TimsVisExtractionFilter),
                                                 POINTER(TimsVisHeatmapSizes)]
    tdf_sdk.tims_vis_calculate_async.restype = c_uint64

    tdf_sdk.tims_vis_cancel.argtypes = [c_uint64]
    tdf_sdk.tims_vis_cancel.restype = c_uint64

    tdf_sdk.tims_vis_close.argtypes = [c_uint64]
    tdf_sdk.tims_vis_close.restype = None

    tdf_sdk.tims_vis_get_chromatogram_line_plot.argtypes = [c_uint64,
                                                            c_int32,
                                                            c_int32,
                                                            c_uint32,
                                                            POINTER(TimsVisLine),
                                                            c_uint32]
    tdf_sdk.tims_vis_get_chromatogram_line_plot.restype = c_uint32

    tdf_sdk.tims_vis_get_image_mob_mz.argtypes = [c_uint64,
                                                  c_uint32,
                                                  POINTER(c_float),
                                                  c_uint32]
    tdf_sdk.tims_vis_get_image_mob_mz.restype = c_uint32

    tdf_sdk.tims_vis_get_image_rt_mob.argtypes = [c_uint64,
                                                  c_uint32,
                                                  POINTER(c_float),
                                                  c_uint32]
    tdf_sdk.tims_vis_get_image_rt_mob.restype = c_uint32

    tdf_sdk.tims_vis_get_image_rt_mz.argtypes = [c_uint64,
                                                 c_uint32,
                                                 POINTER(c_float),
                                                 c_uint32]
    tdf_sdk.tims_vis_get_image_rt_mz.restype = c_uint32

    tdf_sdk.tims_vis_get_last_error_string.argtypes = [c_char_p,
                                                       c_uint32]
    tdf_sdk.tims_vis_get_last_error_string.restype = c_uint32

    tdf_sdk.tims_vis_get_mobilogram_line_plot.argtypes = [c_uint64,
                                                          c_int32,
                                                          c_int32,
                                                          c_uint32,
                                                          POINTER(TimsVisLine),
                                                          c_uint32]
    tdf_sdk.tims_vis_get_mobilogram_line_plot.restype = c_uint32

    tdf_sdk.tims_vis_get_spectrum_line_plot.argtypes = [c_uint64,
                                                        c_int32,
                                                        c_int32,
                                                        c_uint32,
                                                        POINTER(TimsVisLine),
                                                        c_uint32]
    tdf_sdk.tims_vis_get_spectrum_line_plot.restype = c_uint32

    tdf_sdk.tims_vis_get_state.argtypes = [c_uint64,
                                           c_uint64,
                                           c_float,
                                           c_bool]
    tdf_sdk.tims_vis_get_state_restype = c_uint32

    tdf_sdk.tims_vis_open.argtypes = [c_char_p,
                                      c_uint32]
    tdf_sdk.tims_vis_open.restype = c_uint64

    tdf_sdk.tims_vis_wait.argtypes = [c_uint64,
                                      c_uint32]
    tdf_sdk.tims_vis_wait.restype = c_uint32

    tdf_sdk.tims_vis_wait_complete.argtypes = [c_uint64]
    tdf_sdk.tims_vis_wait_complete.restype = c_uint32

    tdf_sdk.tims_voltage_to_scannum.argtypes = convfunc_argtypes
    tdf_sdk.tims_voltage_to_scannum.restype = c_uint32

    tdf_sdk.tsf_close.argtypes = [c_uint64]
    tdf_sdk.tsf_close.restype = None

    tdf_sdk.tsf_get_last_error_string.argtypes = [c_char_p,
                                                  c_uint32]
    tdf_sdk.tsf_get_last_error_string.restype = c_uint32

    tdf_sdk.tsf_has_recalibrated_state.argtypes = [c_uint64]
    tdf_sdk.tsf_has_recalibrated_state.restype = c_uint32

    tdf_sdk.tsf_index_to_mz.argtypes = convfunc_argtypes
    tdf_sdk.tsf_index_to_mz.restype = c_uint32

    tdf_sdk.tsf_mz_to_index.argtypes = convfunc_argtypes
    tdf_sdk.tsf_mz_to_index.restype = c_uint32

    tdf_sdk.tsf_open.argtypes = [c_char_p,
                                 c_uint32]
    tdf_sdk.tsf_open.restype = c_uint64

    tdf_sdk.tsf_read_line_spectrum.argtypes = [c_uint64,
                                               c_int64,
                                               POINTER(c_double),
                                               POINTER(c_float),
                                               c_uint32]
    tdf_sdk.tsf_read_line_spectrum.restype = c_uint32

    tdf_sdk.tsf_read_line_spectrum_v2.argtypes = [c_uint64,
                                                  c_int64,
                                                  POINTER(c_double),
                                                  POINTER(c_float),
                                                  c_int32]
    tdf_sdk.tsf_read_line_spectrum_v2.restype = c_int32

    tdf_sdk.tsf_read_line_spectrum_with_width_v2.argtypes = [c_uint64,
                                                             c_int64,
                                                             POINTER(c_double),
                                                             POINTER(c_float),
                                                             POINTER(c_float),
                                                             c_int32]
    tdf_sdk.tsf_read_line_spectrum_with_width_v2.restype = c_int32

    tdf_sdk.tsf_read_profile_spectrum.argtypes = [c_uint64,
                                                  c_int64,
                                                  POINTER(c_uint32),
                                                  c_uint32]
    tdf_sdk.tsf_read_profile_spectrum.restype = c_uint32

    tdf_sdk.tsf_read_profile_spectrum_v2.argtypes = [c_uint64,
                                                     c_int64,
                                                     POINTER(c_uint32),
                                                     c_int32]
    tdf_sdk.tsf_read_profile_spectrum_v2.restype = c_int32

    tdf_sdk.tsf_set_num_threads.argtypes = [c_uint32]
    tdf_sdk.tsf_set_num_threads.restype = None

    return tdf_sdk
