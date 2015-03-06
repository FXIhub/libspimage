LIST(APPEND OBJECTS ${SPIMAGE_SRC})

IF(CUDA_FOUND)
  LIST(APPEND OBJECTS "fft.cu" "phasing.cu" "phasing_kernels.cu" "support_update_cuda.cu" "image_filter_cuda.cu")
ENDIF(CUDA_FOUND)

IF(BUILD_STATIC_LIB)
	ADD_LIBRARY(spimage STATIC ${OBJECTS})
ELSE(BUILD_STATIC_LIB)
IF(CUDA_FOUND)
IF(CMAKE_COMPILER_IS_GNUCC)
  CUDA_ADD_LIBRARY(spimage SHARED ${OBJECTS} OPTIONS --compiler-options -Wall)
ELSE(CMAKE_COMPILER_IS_GNUCC)
  CUDA_ADD_LIBRARY(spimage SHARED ${OBJECTS})
ENDIF(CMAKE_COMPILER_IS_GNUCC)
  CUDA_ADD_CUFFT_TO_TARGET(spimage)
ELSE(CUDA_FOUND)
  ADD_LIBRARY(spimage SHARED ${OBJECTS})
ENDIF(CUDA_FOUND)
ENDIF(BUILD_STATIC_LIB)



TARGET_LINK_LIBRARIES(spimage ${SPIMAGE_LIBRARIES})


set_target_properties(
 spimage
 PROPERTIES SOVERSION 1
 VERSION 1
 INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib${LIBRARY_SUFFIX}
)



INSTALL(TARGETS spimage
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIBRARY_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIBRARY_SUFFIX})


# Find the python interpreter without calling the findpackage 
# FIND_PROGRAM(PYTHON_EXECUTABLE
#   NAMES python python2.6 python2.5 python2.4 python2.3 python2.2 python2.1 python2.0 python1.6 python1.5
#   PATHS
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.6\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.5\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.4\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.3\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.2\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.1\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.0\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\1.6\\InstallPath]
#   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\1.5\\InstallPath]
#   )

# handle the QUIETLY and REQUIRED arguments and set PYTHONINTERP_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PythonInterp DEFAULT_MSG PYTHON_EXECUTABLE)

MARK_AS_ADVANCED(PYTHON_EXECUTABLE)
if(PYTHONLIBS_FOUND AND SWIG_FOUND AND PYTHONINTERP_FOUND AND PYTHON_NUMPY_FOUND)
#MESSAGE("numpy found compiling python wrappers")
  FILE(GLOB spimage_headers "${CMAKE_CURRENT_SOURCE_DIR}/../include/spimage/*.h")

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c COMMAND ${SWIG_EXECUTABLE} -python -outcurrentdir  -o ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c ${CMAKE_CURRENT_SOURCE_DIR}/spimage_pybackend.i DEPENDS ${spimage_headers} ${CMAKE_CURRENT_SOURCE_DIR}/spimage_pybackend.i)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR})

IF(CUDA_FOUND)
 CUDA_ADD_LIBRARY(_spimage_pybackend SHARED ${OBJECTS} ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c)
 CUDA_ADD_CUFFT_TO_TARGET(_spimage_pybackend)
ELSE(CUDA_FOUND)
 ADD_LIBRARY(_spimage_pybackend SHARED ${OBJECTS} ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c)
ENDIF(CUDA_FOUND)
 TARGET_LINK_LIBRARIES(_spimage_pybackend ${SPIMAGE_LIBRARIES} ${PYTHON_LIBRARIES})
  set_target_properties(
    _spimage_pybackend
    PROPERTIES SOVERSION 1
    VERSION 1
    PREFIX ""
    SUFFIX ".so"
    INSTALL_NAME_DIR ${PYTHON_INSTDIR}
    )
  INSTALL(TARGETS _spimage_pybackend
    RUNTIME DESTINATION ${PYTHON_INSTDIR}
    LIBRARY DESTINATION ${PYTHON_INSTDIR}
    ARCHIVE DESTINATION ${PYTHON_INSTDIR}
    )
  INSTALL(FILES spimage.py ${CMAKE_CURRENT_BINARY_DIR}/spimage_pybackend.py _spimage_conventions.py _spimage_cxi.py _spimage_resample.py _spimage_reconstructor.py _spimage_prtf.py _spimage_io.py _spimage_find_center.py _spimage_sphere_model.py _spimage_spheroid_model.py _spimage_pixelmask.py _spimage_qmap.py _spimage_utils.py _spimage_misc.py _spimage_patterson.py _spimage_material.py elements.dat leastsqbound-scipy/leastsqbound.py DESTINATION  ${PYTHON_INSTDIR})
endif()
