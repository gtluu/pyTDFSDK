import ctypes
import enum


MSMS_SPECTRUM_FUNCTOR = ctypes.CFUNCTYPE(None,
                                         ctypes.c_int64,
                                         ctypes.c_uint32,
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_float))

MSMS_PROFILE_SPECTRUM_FUNCTOR = ctypes.CFUNCTYPE(None,
                                                 ctypes.c_int64,
                                                 ctypes.c_uint32,
                                                 ctypes.POINTER(ctypes.c_int32))
MSMS_SPECTRUM_FUNCTION = ctypes.CFUNCTYPE(None,
                                          ctypes.c_int64,
                                          ctypes.c_uint32,
                                          ctypes.POINTER(ctypes.c_double),
                                          ctypes.POINTER(ctypes.c_float),
                                          ctypes.POINTER(ctypes.c_void_p))

MSMS_PROFILE_SPECTRUM_FUNCTION = ctypes.CFUNCTYPE(None,
                                                  ctypes.c_int64,
                                                  ctypes.c_uint32,
                                                  ctypes.POINTER(ctypes.c_int32),
                                                  ctypes.POINTER(ctypes.c_void_p))


class ChromatogramJob(ctypes.Structure):
    _fields_ = [("id", ctypes.c_int64),
                ("time_begin", ctypes.c_double),
                ("time_end", ctypes.c_double),
                ("mz_min", ctypes.c_double),
                ("mz_max", ctypes.c_double),
                ("ook0_min", ctypes.c_double),
                ("ook0_max", ctypes.c_double)]


CHROMATOGRAM_JOB_GENERATOR = ctypes.CFUNCTYPE(ctypes.c_uint32,
                                              ctypes.POINTER(ChromatogramJob),
                                              ctypes.c_void_p)

CHROMATOGRAM_TRACE_SINK = ctypes.CFUNCTYPE(ctypes.c_uint32,
                                           ctypes.c_int64,
                                           ctypes.c_uint32,
                                           ctypes.POINTER(ctypes.c_int64),
                                           ctypes.POINTER(ctypes.c_uint64),
                                           ctypes.c_void_p)


class TimsVisExtractionFilter(ctypes.Structure):
    _fields_ = [("rtSecondsLower", ctypes.c_double),
                ("rtSecondsUpper", ctypes.c_double),
                ("ook0Lower", ctypes.c_double),
                ("ook0Upper", ctypes.c_double),
                ("mzLower", ctypes.c_double),
                ("mzUpper", ctypes.c_double)]


class TimsVisHeatmapSizes(ctypes.Structure):
    _fields_ = [("widthRtMob", ctypes.c_int32),
                ("heightRtMob", ctypes.c_int32),
                ("widthRtMz", ctypes.c_int32),
                ("heightRtMz", ctypes.c_int32),
                ("widthMobMz", ctypes.c_int32),
                ("heightMobMz", ctypes.c_int32)]


class TimsVisTransformation(enum.Enum):
    NONE = 0
    SQRT = 1
    LIMITED_LOG = 2


class TimsVisLine(ctypes.Structure):
    _fields_ = [("x0", ctypes.c_int32),
                ("y0", ctypes.c_int32),
                ("x1", ctypes.c_int32),
                ("y1", ctypes.c_int32)]


class PressureCompensationStrategy(enum.Enum):
    NoPressureCompensation = 0
    AnalyisGlobalPressureCompensation = 1
    PerFramePressureCompensation = 2
