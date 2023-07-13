import ctypes


def throw_last_tims_data_error(dll_handle):
    len = dll_handle.tims_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(len)
    dll_handle.tims_get_last_error_string(buf, len)
    raise RuntimeError(buf.value)


def throw_last_tsf_data_error (dll_handle):
    len = dll_handle.tsf_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(len)
    dll_handle.tsf_get_last_error_string(buf, len)
    raise RuntimeError(buf.value)
