import os
import platform
import ctypes

if platform.system() == 'Windows':
    TDF_SDK_FILE_NAME = os.path.join(os.path.dirname(__file__), 'tdfsdk2210', 'win64, timsdata.dll')
elif platform.system() == 'Linux':
    TDF_SDK_FILE_NAME = os.path.join(os.path.dirname(__file__), 'tdfsdk2210', 'linux64', 'libtimsdata.so')


def init_tdf_sdk_api(bruker_api_file_name: str=TDF_SDK_FILE_NAME):
    tdf_sdk = ctypes.cdll.LoadLibrary(os.path.realpath(bruker_api_file_name))

    '''
    format
    
    tdf_sdk.[function name].argtypes = [argtypes]
    tdf_sdk.[function name].restype = [return type]
    '''

    # Function 1
    tdf_sdk.tims_ccs_to_oneoverk0_for_mz.argtypes = None
    tdf_sdk.tims_ccs_to_oneoverk0_for_mz.restype = None

    # Function 2
    tdf_sdk.tims_close.argtypes = None
    tdf_sdk.tims_close.restype = None

    # Function 3
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_ext.argtypes = None
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_ext.restype = None

    # Function 4
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_v2.argtypes = None
    tdf_sdk.tims_extract_centroided_spectrum_for_frame_v2.restype = None

    # Function 5
    tdf_sdk.tims_extract_chromatograms.argtypes = None
    tdf_sdk.tims_extract_chromatograms.restype = None

    # Function 6
    tdf_sdk.tims_extract_profile_for_frame.argtypes = None
    tdf_sdk.tims_extract_profile_for_frame.restype = None

    # Function 7
    tdf_sdk.tims_get_last_error_string.argtypes = None
    tdf_sdk.tims_get_last_error_string.restype = None

    # Function 8
    tdf_sdk.tims_has_recalibrated_state.argtypes = None
    tdf_sdk.tims_has_recalibrated_state.restype = None

    # Function 9
    tdf_sdk.tims_index_to_mz.argtypes = None
    tdf_sdk.tims_index_to_mz.restype = None

    # Function 10
    tdf_sdk.tims_mz_to_index.argtypes = None
    tdf_sdk.tims_mz_to_index.restype = None

    # Function 11
    tdf_sdk.tims_oneoverk0_to_ccs_for_mz.argtypes = None
    tdf_sdk.tims_oneoverk0_to_ccs_for_mz.restype = None

    # Function 12
    tdf_sdk.tims_oneoverk0_to_scannum.argtypes = None
    tdf_sdk.tims_oneoverk0_to_scannum.restype = None

    # Function 13
    tdf_sdk.tims_open.argtypes = None
    tdf_sdk.tims_open.restype = None

    # Function 14
    tdf_sdk.tims_open_v2.argtypes = None
    tdf_sdk.tims_open_v2.restype = None

    # Funtion 15
    tdf_sdk.tims_read_pasef_msms.argtypes = None
    tdf_sdk.tims_read_pasef_msms.restype = None

    # Function 16
    tdf_sdk.tims_read_pasef_msms_for_frame.argtypes = None
    tdf_sdk.tims_read_pasef_msms_for_frame.restype = None

    # Function 17
    tdf_sdk.tims_read_pasef_msms_for_frame_v2.argtypes = None
    tdf_sdk.tims_read_pasef_msms_for_frame_v2.restype = None

    # Function 18
    tdf_sdk.tims_read_pasef_msms_v2.argtypes = None
    tdf_sdk.tims_read_pasef_msms_v2.restype = None

    # Function 19
    tdf_sdk.tims_read_pasef_profile_msms.argtypes = None
    tdf_sdk.tims_read_pasef_profile_msms.restype = None

    # Function 20
    tdf_sdk.tims_read_pasef_profile_msms_for_frame.argtypes = None
    tdf_sdk.tims_read_pasef_profile_msms_for_frame.restype = None

    # Function 21
    tdf_sdk.tims_read_pasef_profile_msms_for_frame_v2.argtypes = None
    tdf_sdk.tims_read_pasef_profile_msms_for_frame_v2.restype = None

    # Function 22
    tdf_sdk.tims_read_pasef_profile_msms_v2.argtypes = None
    tdf_sdk.tims_read_pasef_profile_msms_v2.restype = None

    # Function 23
    tdf_sdk.tims_read_scans_internal.argtypes = None
    tdf_sdk.tims_read_scans_internal.restype = None

    # Function 24
    tdf_sdk.tims_read_scans_v2.argtypes = None
    tdf_sdk.tims_read_scans_v2.restype = None

    # Function 25
    tdf_sdk.tims_scannum_to_oneoverk0.argtypes = None
    tdf_sdk.tims_scannum_to_oneoverk0.restype = None

    # Function 26
    tdf_sdk.tims_scannum_to_voltage.argtypes = None
    tdf_sdk.tims_scannum_to_voltage.restype = None

    # Function 27
    tdf_sdk.tims_set_num_threads.argtypes = None
    tdf_sdk.tims_set_num_threads.restype = None

    # Function 28
    tdf_sdk.tims_vis_calculate_async.argtypes = None
    tdf_sdk.tims_vis_calculate_async.restype = None

    # Function 29
    tdf_sdk.tims_vis_cancel.argtypes = None
    tdf_sdk.tims_vis_cancel.restype = None

    # FUnction 30
    tdf_sdk.tims_vis_close.argtypes = None
    tdf_sdk.tims_vis_close.restype = None

    # Function 31
    tdf_sdk.tims_vis_get_chromatogram_line_plot.argtypes = None
    tdf_sdk.tims_vis_get_chromatogram_line_plot.restype = None

    # Function 32
    tdf_sdk.tims_vis_get_image_mob_mz.argtypes = None
    tdf_sdk.tims_vis_get_image_mob_mz.restype = None

    # Function 33
    tdf_sdk.tims_vis_get_image_rt_mob.argtypes = None
    tdf_sdk.tims_vis_get_image_rt_mob.restype = None

    # Function 34
    tdf_sdk.tims_vis_get_image_rt_mz.argtypes = None
    tdf_sdk.tims_vis_get_image_rt_mz.restype = None

    # Function 35
    tdf_sdk.tims_vis_get_last_error_string.argtypes = None
    tdf_sdk.tims_vis_get_last_error_string.restype = None

    # Function 36
    tdf_sdk.tims_vis_get_mobilogram_line_plot.argtypes = None
    tdf_sdk.tims_vis_get_mobilogram_line_plot.restype = None

    # Function 37
    tdf_sdk.tims_vis_get_spectrum_line_plot.argtypes = None
    tdf_sdk.tims_vis_get_spectrum_line_plot.restype = None

    # Function 38
    tdf_sdk.tims_vis_get_state.argtypes = None
    tdf_sdk.tims_vis_get_state_restype = None

    # Function 39
    tdf_sdk.tims_vis_open.argtypes = None
    tdf_sdk.tims_vis_open.restype = None

    # Function 40
    tdf_sdk.tims_vis_wait.argtypes = None
    tdf_sdk.tims_vis_wait.restype = None

    # Function 41
    tdf_sdk.tims_vis_wait_complete.argtypes = None
    tdf_sdk.tims_vis_wait_complete.restype = None

    # Function 42
    tdf_sdk.tims_voltage_to_scannum.argtypes = None
    tdf_sdk.tims_voltage_to_scannum.restype = None

    # Function 43
    tdf_sdk.tsf_close.argtypes = None
    tdf_sdk.tsf_close.restype = None

    # Function 44
    tdf_sdk.tsf_get_last_error_string.argtypes = None
    tdf_sdk.tsf_get_last_error_string.restype = None

    # Function 45
    tdf_sdk.tsf_has_recalibrated_state.argtypes = None
    tdf_sdk.tsf_has_recalibrated_state.restype = None

    # Function 46
    tdf_sdk.tsf_index_to_mz.argtypes = None
    tdf_sdk.tsf_index_to_mz.restype = None

    # Function 47
    tdf_sdk.tsf_mz_to_index.argtypes = None
    tdf_sdk.tsf_mz_to_index.restype = None

    # Function 48
    tdf_sdk.tsf_open.argtypes = None
    tdf_sdk.tsf_open.restype = None

    # Function 49
    tdf_sdk.tsf_read_line_spectrum.argtypes = None
    tdf_sdk.tsf_read_line_spectrum.restype = None

    # Function 50
    tdf_sdk.tsf_read_line_spectrum.argtypes = None
    tdf_sdk.tsf_read_line_spectrum.restype = None

    # Function 51
    tdf_sdk.tsf_read_line_spectrum_with_width.argtypes = None
    tdf_sdk.tsf_read_line_spectrum_with_width.restype = None

    # Function 52
    tdf_sdk.tsf_read_line_spectrum_with_width_v2.argtypes = None
    tdf_sdk.tsf_read_line_spectrum_with_width_v2.restype = None

    # Function 53
    tdf_sdk.tsf_read_profile_spectrum.argtypes = None
    tdf_sdk.tsf_read_profile_spectrum.restype = None

    # Function 54
    tdf_sdk.tsf_read_profile_spectrum_v2.argtypes = None
    tdf_sdk.tsf_read_profile_spectrum_v2.restype = None

    # Function 55
    tdf_sdk.tsf_set_num_threads.argtypes = None
    tdf_sdk.tsf_set_num_threads.restype = None
    