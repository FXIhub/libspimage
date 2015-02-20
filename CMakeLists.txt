cmake_minimum_required(VERSION 2.6)

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # CMake 3.0
endif()

if(POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW) # CMake 3.0
endif()

PROJECT (spimage)

INCLUDE(CTest)
ENABLE_TESTING()
IF(BUILD_TESTING)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)


SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Make sure that while building this library the include files are not 
# built with export prefix
ADD_DEFINITIONS(-D SPIMAGE_NO_DLL)

IF(WIN32)
#	ADD_DEFINITIONS(-D_HDF5USEDLL_)
	ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
ENDIF(WIN32)

IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # Necessary due to CUDA compatibility problems
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")  
  # The CUDA module thinks cc is gcc and passes dumpspecs.
  # This is a workaround
  if(CUDA_HOST_COMPILER STREQUAL "/usr/bin/cc")
    SET(CUDA_HOST_COMPILER "/usr/bin/clang")
  endif()
endif()

	

SET(BUILD_STATIC_LIB CACHE BOOL "If ON output static library. Otherwise output shared libraries.")
SET(INCLUDE_DEPENDENCIES CACHE BOOL "[EXPERIMENTAL] If ON embed dependent libraries.")
SET(DMALLOC_USE OFF CACHE BOOL "If ON link to dmalloc library if possible.")
SET(DOUBLE_PRECISION OFF CACHE BOOL "If ON use double precision. Otherwise use single precision")
SET(SP_MEM_DEBUG OFF CACHE BOOL "If ON use memory debugging code.")
SET(USE_CUDA ON CACHE BOOL "If ON try to use CUDA.")
SET(PYTHON_WRAPPERS ON CACHE BOOL "If ON try to build python wrappers.")

find_package(TIFF)
find_package(FFTW3 REQUIRED)
find_package(PNG)
find_package(HDF5 REQUIRED)
find_package(GSL)

if(PYTHON_WRAPPERS)
  find_package(PythonLibs)
  find_package(PythonInterp)
  find_package(Numpy)
  find_package(SWIG)
  if(PYTHONLIBS_FOUND AND SWIG_FOUND AND PYTHONINTERP_FOUND AND PYTHON_NUMPY_FOUND)
    # Taken from plplot. Many thanks!
    # N.B. This is a nice way to obtain all sorts of python information
    # using distutils.
    execute_process(
	 COMMAND
  	 ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_python_lib()"
  	 OUTPUT_VARIABLE _PYTHON_INSTDIR
  	 OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    SET(PYTHON_INSTDIR ${_PYTHON_INSTDIR} CACHE PATH "Installation directory for python module.")
  endif()
endif(PYTHON_WRAPPERS)
IF(USE_CUDA)
  IF(CMAKE_VERSION)
    # CMAKE VERSION was first defined in 2.6.3. Sorry for the hack
    IF(CMAKE_VERSION VERSION_LESS 2.8)
      SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/cuda")
    ELSE(CMAKE_VERSION VERSION_LESS 2.8)
      MESSAGE(STATUS "Using builtin CUDA package")
    ENDIF(CMAKE_VERSION VERSION_LESS 2.8)
    find_package(CUDA QUIET)
  ELSE(CMAKE_VERSION)
    message(STATUS "cmake 2.6.3 or higher required to use CUDA")
  ENDIF(CMAKE_VERSION)

  IF(CUDA_FOUND)
    mark_as_advanced(CUDA_SDK_ROOT_DIR)
    mark_as_advanced(CUDA_BUILD_EMULATION)
    mark_as_advanced(CUDA_BUILD_CUBIN)
    mark_as_advanced(CUDA_VERBOSE_BUILD)
    ADD_DEFINITIONS(-D_USE_CUDA)
    # Set a reasonable minimum arch
    if(NOT CUDA_NVCC_FLAGS)
      # This ensure that the code also runs on cards with computing cap. as low as 1.1
      SET(CUDA_NVCC_FLAGS "-gencode=arch=compute_11,code=sm_11 -gencode=arch=compute_20,code=sm_20")
      #SET(CUDA_NVCC_FLAGS "-arch=sm_20")
    endif()
    message(STATUS "CUDA found. Using CUDA.")
  ELSE(CUDA_FOUND)
    message(STATUS "CUDA not found. Not using CUDA.")
  ENDIF(CUDA_FOUND)
ELSE(USE_CUDA)
  UNSET(CUDA_FOUND)
ENDIF(USE_CUDA)

#   Offer the user the choice of overriding the installation directories

get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

if (${LIB64} STREQUAL "TRUE")
    set(LIBRARY_SUFFIX  64)
else()
    set(LIBRARY_SUFFIX  "")
endif()



IF(DOUBLE_PRECISION)
	ADD_DEFINITIONS(-D_SP_DOUBLE_PRECISION)
ELSE(DOUBLE_PRECISION)	
	ADD_DEFINITIONS(-D_SP_SINGLE_PRECISION)
ENDIF(DOUBLE_PRECISION)

IF(SP_MEM_DEBUG)
	ADD_DEFINITIONS(-D_SP_MEM_DEBUG)
ENDIF(SP_MEM_DEBUG)
	

IF(DMALLOC_USE)
  FIND_LIBRARY (DMALLOC_LIBRARY dmalloc /usr/lib /sw/lib)
  FIND_PATH (DMALLOC_INCLUDE_DIR dmalloc.h PATHS /usr/include /sw/include)	
  IF(DMALLOC_LIBRARY)
    ADD_DEFINITIONS(-D_USE_DMALLOC)
    SET(LINK_TO_DMALLOC 1)
  ENDIF(DMALLOC_LIBRARY)
ENDIF(DMALLOC_USE)


# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

IF(APPLE)
  # There's a problem with install_name_tool which cannot
  # handle duplicate rpaths, so we'll take the calculated risk
  # of not including rpaths in the build on macosx
  SET(CMAKE_SKIP_BUILD_RPATH  TRUE)
ELSE(APPLE)
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
ENDIF(APPLE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib${LIBRARY_SUFFIX}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIBRARY_SUFFIX}")
ENDIF("${isSystemDir}" STREQUAL "-1")


configure_file(${PROJECT_SOURCE_DIR}/include/config.h.in 
  ${PROJECT_BINARY_DIR}/config.h)
include_directories("${PROJECT_BINARY_DIR}")


INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR} ${PNG_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${GSL_INCLUDE_DIR})
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_SOURCE_DIR}/include/")
IF(CUDA_FOUND)
INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
ENDIF(CUDA_FOUND)

LIST(APPEND TESTS_LIBRARIES  ${TIFF_LIBRARIES} ${FFTW3_LIBRARIES} ${PNG_LIBRARIES} ${HDF5_LIBRARIES}  ${GSL_LIBRARIES})
LIST(APPEND SPIMAGE_LIBRARIES  ${TIFF_LIBRARIES} ${FFTW3_LIBRARIES} ${PNG_LIBRARIES} ${HDF5_LIBRARIES})
	    
IF(LINK_TO_DMALLOC)
LIST(APPEND SPIMAGE_LIBRARIES ${DMALLOC_LIBRARY})
LIST(APPEND TESTS_LIBRARIES ${DMALLOC_LIBRARY})
ENDIF(LINK_TO_DMALLOC)

LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/image_util.c" "${CMAKE_SOURCE_DIR}/src/fft.c" "${CMAKE_SOURCE_DIR}/src/linear_alg.c")
LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/mem_util.c" "${CMAKE_SOURCE_DIR}/src/image_sphere.c" "${CMAKE_SOURCE_DIR}/src/sperror.c")
LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/gaussianinv.c" "${CMAKE_SOURCE_DIR}/src/image_noise.c" "${CMAKE_SOURCE_DIR}/src/hashtable.c")
LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/interpolation_kernels.c" "${CMAKE_SOURCE_DIR}/src/time_util.c")
LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/list.c" "${CMAKE_SOURCE_DIR}/src/prtf.c" "${CMAKE_SOURCE_DIR}/src/phasing.c")
LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/colormap.c" "${CMAKE_SOURCE_DIR}/src/cuda_util.c" "${CMAKE_SOURCE_DIR}/src/support_update.c")
LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/image_io.c" "${CMAKE_SOURCE_DIR}/src/image_filter.c")
LIST(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/find_center.c")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(include)

if(GSL_FOUND)
ADD_SUBDIRECTORY(tests)
endif()
