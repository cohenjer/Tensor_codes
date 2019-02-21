import ctypes as ct
from ctypes import c_int, c_double, c_void_p, byref, POINTER
import numpy as np
from numpy import int32, float64

libfg = ct.CDLL('libfgnsr.so')

# Shorthands for array types we use for argument checking
dbl_arr = np.ctypeslib.ndpointer(dtype=float64, ndim=1, flags='C_CONTIGUOUS')
int_arr = np.ctypeslib.ndpointer(dtype=int32,   ndim=1, flags='C_CONTIGUOUS')

def handle_status(status):
    if status.value != 0:
        raise Error("Implement a better error handling today!")


libfg.FGNSRproj_spinit.argtypes = [
    c_int,             # lenvec
    dbl_arr,           # weights
    c_int,             # verbose
    POINTER(c_void_p), # data pointer (by pointer)
] 
libfg.FGNSRproj_spinit.restype = c_int
def FGNSRproj_spinit(lenvec, weights, verbose):
    data_p = c_void_p()
    status = libfg.FGNSRproj_spinit(lenvec, weights, verbose, byref(data_p))
    handle_status(status)
    return data_p

libfg.FGNSRproj_spfree.argtypes = [
    POINTER(c_void_p), # data pointer (by pointer)
]
def FGNSRproj_spfree(data_ptr):
    libfg.FGNSRproj_spfree(byref(data_ptr))

libfg.FGNSRproj_project_spmatrix.argtypes = [
    POINTER(c_void_p), # data pointer
    c_int,             # numvec
    int_arr,           # begin
    int_arr,           # index
    dbl_arr,           # value
    c_double,          # ub
    int_arr,           # diag_idx
]
libfg.FGNSRproj_project_spmatrix.restype = c_int

def FGNSRproj_project_spmatrix(data_p, numvec, begin, index, value, ub,
        diag_idx):
    pass

