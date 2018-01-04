[![Build Status](https://travis-ci.org/FXIhub/libspimage.svg?branch=master)](https://travis-ci.org/FXIhub/libspimage)

LibSPImage - Image Processing Library
=====================================

LibSPImage is a library for image processing used by Hawk.

<br>
Authors:<br>
Filipe Maia     ||  <filipe.c.maia@gmail.com><br>
Tomas Ekeberg   ||  <ekeberg@xray.bmc.uu.se><br>
Max Hantke      ||  <max.hantke@icm.uu.se><br>
Benedikt Daurer ||  <benedikt@xray.bmc.uu.se><br>
Jonas Sellberg  ||  <sellberg@xray.bmc.uu.se><br>
<br>
LibSPImage is currently v1.00 and is undergoing testing.

-------------------------------------------------------------------------------


Dependencies
------------

LibSPImage depends on the following packages:

* cmake
* libpng
* libtiff
* zlib
* szip
* HDF5 (v1.8)
* FFTW (v3.x)
* GSL

The packages above are necessary for the compilation to complete. Additionally, LibSPImage also uses the following packages for the Python wrapper, including the Reconstruction class and the prtf function:

* pcre
* SWIG
* python (v2.7 | v3.4 | v3.5 | v3.6)
* numpy
* h5py

We recommend you to install all these prerequisites before continuing any further. If you're unfamiliar with installing packages, we recommend you to use Homebrew, which can be found at http://brew.sh/

If you want to be able to run LibSPImage on your graphics card, you need CUDA, which can be downloaded at https://developer.nvidia.com/cuda-downloads


Installation
------------

Once all the dependencies are installed, building and installing LibSPImage should be completed in a few simple steps:

- Clone the repository to your local computer:

    `git clone https://github.com/filipemaia/libspimage.git`

If you don't have git installed, you can follow the link and click `Download ZIP`, but we recommend that you install git (available through Homebrew).

- Create and go into a build directory:

    `cd libspimage`

    `mkdir build`

    `cd build`

- Run ccmake and point it to the base directory:

    `ccmake ..`

You will see something like:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[EMPTY CACHE]

EMPTY CACHE:

Press [enter] to edit option                              CMake Version 3.0.2
Press [c] to configure
Press [h] for help               Press [q] to quit without generating
Press [t] to toggle advanced mode (Currently Off)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Press `c` to configure

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BUILD_STATIC_LIB                *OFF                                                                                                                
BUILD_TESTING                   *ON                                                                                                                 
CMAKE_BUILD_TYPE                *                                                                                                                   
CMAKE_INSTALL_PREFIX            */usr/local                                                                                                         
CMAKE_OSX_ARCHITECTURES         *                                                                                                                   
CMAKE_OSX_DEPLOYMENT_TARGET     *                                                                                                                   
CMAKE_OSX_SYSROOT               */Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk                 
CUDA_HOST_COMPILER              */usr/bin/cc                                                                                                        
CUDA_SEPARABLE_COMPILATION      *OFF                                                                                                                
CUDA_TOOLKIT_ROOT_DIR           */usr/local/cuda                                                                                                    
DMALLOC_USE                     *OFF                                                                                                                
DOUBLE_PRECISION                *OFF                                                                                                                
GSL_CONFIG                      */usr/local/bin/gsl-config                                                                                          
GSL_CONFIG_PREFER_PATH          */bin;;/bin;                                                                                                        
GSL_EXE_LINKER_FLAGS            *-Wl,-rpath,/usr/local/Cellar/gsl/1.16/lib                                                                          
HDF5_INCLUDE_DIR                */usr/local/include                                                                                                 
HDF5_LIBRARY                    */usr/local/lib/libhdf5.dylib                                                                                       
INCLUDE_DEPENDENCIES            *OFF                                                                                                                
PYTHON_INSTDIR                  */Library/Python/2.7/site-packages                                                                                  
PYTHON_WRAPPERS                 *ON                                                                                                                 
SP_MEM_DEBUG                    *OFF                                                                                                                
SWIG_EXECUTABLE                 */usr/local/bin/swig                                                                                                
USE_CUDA                        *ON                                                                                                                 

Press [enter] to edit option                              CMake Version 3.0.2
Press [c] to configure
Press [h] for help               Press [q] to quit without generating
Press [t] to toggle advanced mode (Currently Off)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure that all the dependencies have been found if they're installed in non-standard paths. You can also specify `CMAKE_INSTALL_PREFIX` to set the install directory. This is done by highlighting `CMAKE_INSTALL_PREFIX`, pressing enter to modify the path, then press enter again to set it. You may press `t` to view all available options.

- If you made any changes, press `c` again to re-configure.

- If everything went well you should see a screen just like the one above and be able to press `g` to generate the Makefiles and exit.

- Now run:

    `make`

This will build things and place the result in the build directory. If the build is slow, you can specify `make -j 2` to use 2 threads.

- If you want to install, just run:

    `make install`

This will automatically install the API in your python distribution and put the library and its header files in your library path. You may have to obtain administrative permissions by typing `sudo make install` for the install to successfully complete.


Help
----

There are two sources of documentation for LibSPImage:

* This README
* The LibSPImage API is documented in-place in the code


Contribute
----------

If you would like to add content to LibSPImage or suggest an enhancement, please do! The usual GitHub features (issues & pull requests) should suffice for these purposes.
