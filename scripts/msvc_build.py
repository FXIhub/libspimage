#!/usr/bin/env python

env_map = os.environ;

p = subprocess.Popen(['cmake',
                  '-G','Visual Studio 12 Win64',
                  '-DTIFF_LIBRARY:FILEPATH=c:/Program Files (x86)/GnuWin32/lib/libtiff.lib',
                  '-DTIFF_INCLUDE_DIR:PATH=c:/Program Files (x86)/GnuWin32/include',
                  '-DJPEG_LIBRARY:FILEPATH=c:/Program Files (x86)/GnuWin32/bin/jpeg.lib',
                  '-DJPEG_INCLUDE_DIR:PATH=c:/Program Files (x86)/GnuWin32/include',
                  '-DPNG_LIBRARY:FILEPATH=c:/Program Files (x86)/GnuWin32/bin/libpng.lib',
                  '-DPNG_PNG_INCLUDE_DIR:PATH=c:/Program Files (x86)/GnuWin32/include',
                  '-DZLIB_LIBRARY:FILEPATH=c:/Program Files (x86)/GnuWin32/bin/zlib.lib',
                  '-DZLIB_INCLUDE_DIR:PATH=c:/Program Files (x86)/GnuWin32/include',
                  '-DFFTW3_FFTWF_LIBRARY:FILEPATH=c:/fftw3/libfftw3f-3.lib',
                  '-DFFTW3_FFTW_LIBRARY:FILEPATH=c:/fftw3/libfftw3-3.lib',
                  '-DFFTW3_FFTWL_LIBRARY:FILEPATH=c:/fftw3/libfftw3l-3.lib',
                  '-DFFTW3_INCLUDE_DIR:PATH=c:/fftw3',
                  '-DHDF5_LIBRARY:FILEPATH=c:/Program Files (x86)/HDF_Group/HDF5/1.8.15/lib/hdf5.lib',
                  '-DHDF5_INCLUDE_DIR:PATH=c:/Program Files (x86)/HDF_Group/HDF5/1.8.15/include',
                  '-DGSL_GSL_LIBRARY:FILEPATH=c:/Program Files (x86)/GnuWin32/lib/libgsl.lib',
                  '-DGSL_GSLCBLAS_LIBRARY:FILEPATH=c:/Program Files (x86)/GnuWin32/lib/libgslcblas.lib',
                  '-DGSL_INCLUDE_DIR:PATH=c:/Program Files (x86)/GnuWin32/include',
                  '-DCMAKE_VERBOSE_MAKEFILE:BOOL=FALSE',
                  '..'], env=env_map)

