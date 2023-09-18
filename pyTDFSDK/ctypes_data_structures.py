from ctypes import CFUNCTYPE, POINTER, Structure, c_int64, c_uint32, c_double, c_float, c_int32, c_void_p, c_uint64
from enum import Enum

MSMS_SPECTRUM_FUNCTOR = CFUNCTYPE(None,
                                  c_int64,
                                  c_uint32,
                                  POINTER(c_double),
                                  POINTER(c_float))

MSMS_PROFILE_SPECTRUM_FUNCTOR = CFUNCTYPE(None,
                                          c_int64,
                                          c_uint32,
                                          POINTER(c_int32))
MSMS_SPECTRUM_FUNCTION = CFUNCTYPE(None,
                                   c_int64,
                                   c_uint32,
                                   POINTER(c_double),
                                   POINTER(c_float),
                                   POINTER(c_void_p))

MSMS_PROFILE_SPECTRUM_FUNCTION = CFUNCTYPE(None,
                                           c_int64,
                                           c_uint32,
                                           POINTER(c_int32),
                                           POINTER(c_void_p))


class ChromatogramJob(Structure):
    """
    ctypes.Structure to hold chromatogram job information.
    """
    _fields_ = [("id", c_int64),
                ("time_begin", c_double),
                ("time_end", c_double),
                ("mz_min", c_double),
                ("mz_max", c_double),
                ("ook0_min", c_double),
                ("ook0_max", c_double)]


CHROMATOGRAM_JOB_GENERATOR = CFUNCTYPE(c_uint32,
                                       POINTER(ChromatogramJob),
                                       c_void_p)

CHROMATOGRAM_TRACE_SINK = CFUNCTYPE(c_uint32,
                                    c_int64,
                                    c_uint32,
                                    POINTER(c_int64),
                                    POINTER(c_uint64),
                                    c_void_p)


class TimsVisExtractionFilter(Structure):
    """
    ctypes.Structure to hold visualization extraction filter information.
    """
    _fields_ = [("rtSecondsLower", c_double),
                ("rtSecondsUpper", c_double),
                ("ook0Lower", c_double),
                ("ook0Upper", c_double),
                ("mzLower", c_double),
                ("mzUpper", c_double)]


class TimsVisHeatmapSizes(Structure):
    """
    ctypes.Structure to hold heatmap size information.
    """
    _fields_ = [("widthRtMob", c_int32),
                ("heightRtMob", c_int32),
                ("widthRtMz", c_int32),
                ("heightRtMz", c_int32),
                ("widthMobMz", c_int32),
                ("heightMobMz", c_int32)]


class TimsVisTransformation(Enum):
    """
    enum.Enum to determine MSI visualization intensity transformation method.
    """
    NONE = 0
    SQRT = 1
    LIMITED_LOG = 2


class TimsVisLine(Structure):
    """
    ctypes.Structure to hold visualization x-y coordinates.
    """
    _fields_ = [("x0", c_int32),
                ("y0", c_int32),
                ("x1", c_int32),
                ("y1", c_int32)]


class PressureCompensationStrategy(Enum):
    """
    enum.Enum to determine pressure compensation strategy.
    """
    NoPressureCompensation = 0
    AnalyisGlobalPressureCompensation = 1
    PerFramePressureCompensation = 2
