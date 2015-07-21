#!/usr/bin/env python

import os
import subprocess

env_map = os.environ;

p = subprocess.Popen(['cmake',
                  '-G','Visual Studio 12 Win64',
                  '-DTIFF_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libtiff.lib',
                  '-DTIFF_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DPNG_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libpng.lib',
                  '-DPNG_PNG_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DZLIB_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libz.lib',
                  '-DZLIB_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DFFTW3_FFTWF_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libfftw3f.lib',
                  '-DFFTW3_FFTW_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libfftw3.lib',
                  '-DFFTW3_FFTWL_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libfftw3l.lib',
                  '-DFFTW3_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DHDF5_LIBRARY:FILEPATH=c:/Program Files/HDF_Group/HDF5/1.8.15/lib/hdf5.lib',
                  '-DHDF5_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DGSL_GSL_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libgsl.lib',
                  '-DGSL_GSLCBLAS_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libgslcblas.lib',
                  '-DGSL_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DCMAKE_VERBOSE_MAKEFILE:BOOL=FALSE',
                  '..'], env=env_map)

