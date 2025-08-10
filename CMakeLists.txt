
cmake_minimum_required(VERSION 3.18)

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # CMake 3.0
endif()

if(POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW) # CMake 3.0
endif()

project(spimage)

include(CTest)
enable_testing()
if(BUILD_TESTING)
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)
endif()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Make sure that while building this library the include files are not 
# built with export prefix
add_definitions(-D SPIMAGE_NO_DLL)

if(WIN32)
#   add_definitions(-D_HDF5USEDLL_)
add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
# Required for compilign with MSVC as documented in
# USING_HDF5_VS.txt from the hdf group
add_definitions(-DH5_BUILT_AS_DYNAMIC_LIB)
endif()

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
endif()

  

set(BUILD_STATIC_LIB CACHE BOOL "If ON output static library. Otherwise output shared libraries.")
set(INCLUDE_DEPENDENCIES CACHE BOOL "[EXPERIMENTAL] If ON embed dependent libraries.")
set(DMALLOC_USE OFF CACHE BOOL "If ON link to dmalloc library if possible.")
set(DOUBLE_PRECISION OFF CACHE BOOL "If ON use double precision. Otherwise use single precision")
set(SP_MEM_DEBUG OFF CACHE BOOL "If ON use memory debugging code.")
set(USE_CUDA ON CACHE BOOL "If ON try to use CUDA.")
set(PYTHON_WRAPPERS ON CACHE BOOL "If ON try to build python wrappers.")

find_package(TIFF)
find_package(FFTW3 REQUIRED)
find_package(PNG)
find_package(HDF5 REQUIRED)
find_package(GSL)



if(PYTHON_WRAPPERS)
  find_package (Python3 COMPONENTS Interpreter Development NumPy)
  if(Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE PATH "Python interpreter executable")
    set(PYTHON_INCLUDE_PATH ${Python3_INCLUDE_DIRS} CACHE PATH "Python include dirs")
    set(PYTHON_LIBRARIES ${Python3_LIBRARIES} CACHE PATH "Python libraries")
    set(PYTHON_NUMPY_INCLUDE_DIR ${Python3_NumPy_INCLUDE_DIRS} CACHE PATH "Python NumPy include dirs")
    set(PYTHON_NUMPY_FOUND ${Python3_NumPy_FOUND})
    mark_as_advanced(PYTHON_INCLUDE_PATH PYTHON_LIBRARIES PYTHON_NUMPY_INCLUDE_DIR)
  else()
    message(FATAL_ERROR "Python interpreter not found")
  endif()
  find_package(SWIG)
  if(PYTHON_LIBRARIES AND SWIG_FOUND AND PYTHON_EXECUTABLE AND PYTHON_NUMPY_FOUND)
    set(PYTHON_INSTDIR ${Python3_SITELIB} CACHE PATH "Installation directory for python module.")
    mark_as_advanced(PYTHON_INSTDIR)
  endif()
endif()
if(USE_CUDA)
  find_package(CUDAToolkit)
  if(CUDAToolkit_BIN_DIR)
    set(CUDA_FOUND TRUE)
  else()
    set(CUDA_FOUND FALSE)
  endif()
  if(CUDA_FOUND)
    enable_language(CUDA)
    add_definitions(-D_USE_CUDA)
    # Set a reasonable minimum arch
    if(CUDAToolkit_VERSION_MAJOR GREATER 11)
      set(CMAKE_CUDA_ARCHITECTURES "50" CACHE STRING "CUDA architectures to build for")
    elseif(CUDAToolkit_VERSION_MAJOR GREATER 10)
      set(CMAKE_CUDA_ARCHITECTURES "35" CACHE STRING "CUDA architectures to build for")
    else()
      set(CMAKE_CUDA_ARCHITECTURES "30" CACHE STRING "CUDA architectures to build for")
    endif()
    message(STATUS "CUDA found. Using CUDA.")
  else(CUDA_FOUND)
    message(STATUS "CUDA not found. Not using CUDA.")
  endif(CUDA_FOUND)
else(USE_CUDA)
  unset(CUDA_FOUND)
endif(USE_CUDA)

#   Offer the user the choice of overriding the installation directories
include(GNUInstallDirs)
if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

if(DOUBLE_PRECISION)
  add_definitions(-D_SP_DOUBLE_PRECISION)
else(DOUBLE_PRECISION)
  add_definitions(-D_SP_SINGLE_PRECISION)
endif(DOUBLE_PRECISION)

if(SP_MEM_DEBUG)
  add_definitions(-D_SP_MEM_DEBUG)
endif(SP_MEM_DEBUG)

if(DMALLOC_USE)
  find_library(DMALLOC_LIBRARY dmalloc /usr/lib /sw/lib)
  find_path(DMALLOC_INCLUDE_DIR dmalloc.h PATHS /usr/include /sw/include)
  if(DMALLOC_LIBRARY)
    add_definitions(-D_USE_DMALLOC)
    set(LINK_TO_DMALLOC 1)
  endif(DMALLOC_LIBRARY)
endif(DMALLOC_USE)


# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

if(APPLE)
  # There's a problem with install_name_tool which cannot
  # handle duplicate rpaths, so we'll take the calculated risk
  # of not including rpaths in the build on macosx
  set(CMAKE_SKIP_BUILD_RPATH  TRUE)
else(APPLE)
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)
endif(APPLE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${CMAKE_INSTALL_LIBDIR} isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_LIBDIR})
endif("${isSystemDir}" STREQUAL "-1")


configure_file(${PROJECT_SOURCE_DIR}/include/config.h.in 
  ${PROJECT_BINARY_DIR}/config.h)
include_directories("${PROJECT_BINARY_DIR}")


include_directories(${FFTW3_INCLUDE_DIR} ${PNG_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${GSL_INCLUDE_DIR})
include_directories(BEFORE "${CMAKE_SOURCE_DIR}/include/")
if(CUDA_FOUND)
include_directories(${CUDA_INCLUDE_DIRS})
endif()


list(APPEND TESTS_LIBRARIES  ${TIFF_LIBRARIES} ${FFTW3_LIBRARIES} ${PNG_LIBRARIES} ${HDF5_LIBRARIES}  ${GSL_LIBRARIES})
list(APPEND SPIMAGE_LIBRARIES  ${TIFF_LIBRARIES} ${FFTW3_LIBRARIES} ${PNG_LIBRARIES} ${HDF5_LIBRARIES})
    
if(LINK_TO_DMALLOC)
list(APPEND SPIMAGE_LIBRARIES ${DMALLOC_LIBRARY})
list(APPEND TESTS_LIBRARIES ${DMALLOC_LIBRARY})
endif()

list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/image_util.c" "${CMAKE_SOURCE_DIR}/src/fft.c" "${CMAKE_SOURCE_DIR}/src/linear_alg.c")
list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/mem_util.c" "${CMAKE_SOURCE_DIR}/src/image_sphere.c" "${CMAKE_SOURCE_DIR}/src/sperror.c")
list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/gaussianinv.c" "${CMAKE_SOURCE_DIR}/src/image_noise.c" "${CMAKE_SOURCE_DIR}/src/hashtable.c")
list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/interpolation_kernels.c" "${CMAKE_SOURCE_DIR}/src/time_util.c")
list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/list.c" "${CMAKE_SOURCE_DIR}/src/prtf.c" "${CMAKE_SOURCE_DIR}/src/phasing.c")
list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/colormap.c" "${CMAKE_SOURCE_DIR}/src/cuda_util.c" "${CMAKE_SOURCE_DIR}/src/support_update.c")
list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/image_io.c" "${CMAKE_SOURCE_DIR}/src/image_filter.c")
list(APPEND SPIMAGE_SRC "${CMAKE_SOURCE_DIR}/src/find_center.c")

add_subdirectory(src)
add_subdirectory(include)

if(GSL_FOUND)
add_subdirectory(tests)
endif()
