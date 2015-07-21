#!/usr/bin/env python

import os
import subprocess

env_map = os.environ;

p = subprocess.Popen(['cmake',
                  '-G','Visual Studio 12 Win64',
                  '-DTIFF_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libtiff.dll.a',
                  '-DTIFF_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DPNG_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libpng.dll.a',
                  '-DPNG_PNG_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DZLIB_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libz.dll.a',
                  '-DZLIB_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DFFTW3_FFTWF_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libfftw3f.dll.a',
                  '-DFFTW3_FFTW_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libfftw3.dll.a',
                  '-DFFTW3_FFTWL_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libfftw3l.dll.a',
                  '-DFFTW3_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DHDF5_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libhdf5.dll.a',
                  '-DHDF5_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DGSL_GSL_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libgsl.dll.a',
                  '-DGSL_GSLCBLAS_LIBRARY:FILEPATH=c:/msys64/mingw64/lib/libgslcblas.dll.a',
                  '-DGSL_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DCMAKE_VERBOSE_MAKEFILE:BOOL=FALSE',
                  '..'], env=env_map)

