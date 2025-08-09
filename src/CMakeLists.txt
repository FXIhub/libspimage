LIST(APPEND OBJECTS ${SPIMAGE_SRC})



IF(CUDA_FOUND)
  LIST(APPEND OBJECTS "fft.cu" "phasing.cu" "phasing_kernels.cu" "support_update_cuda.cu" "image_filter_cuda.cu")
ENDIF(CUDA_FOUND)

IF(BUILD_STATIC_LIB)
	ADD_LIBRARY(spimage STATIC ${OBJECTS})
ELSE(BUILD_STATIC_LIB)
ADD_LIBRARY(spimage SHARED ${OBJECTS})
IF(CUDA_FOUND)
  LIST(APPEND SPIMAGE_LIBRARIES CUDA::curand CUDA::cufft)
  #set_target_properties(spimage PROPERTIES CUDA_ARCHITECTURES "35")
ENDIF(CUDA_FOUND)
ENDIF(BUILD_STATIC_LIB)

TARGET_LINK_LIBRARIES(spimage PRIVATE ${SPIMAGE_LIBRARIES})

set_target_properties(
 spimage
 PROPERTIES
 SOVERSION 1
 VERSION 1
 INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR}
)

if(WIN32)
  set_target_properties(
    spimage
    PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif(WIN32)


INSTALL(TARGETS spimage
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(SWIG_FOUND AND Python3_Interpreter_FOUND AND PYTHON_NUMPY_FOUND)
  FILE(GLOB spimage_headers "${CMAKE_CURRENT_SOURCE_DIR}/../include/spimage/*.h")

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c COMMAND ${SWIG_EXECUTABLE} -python -outcurrentdir  -o ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c ${CMAKE_CURRENT_SOURCE_DIR}/spimage_pybackend.i DEPENDS ${spimage_headers} ${CMAKE_CURRENT_SOURCE_DIR}/spimage_pybackend.i)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR})

  ADD_LIBRARY(_spimage_pybackend SHARED ${OBJECTS} ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c)

  if(APPLE)
    TARGET_LINK_LIBRARIES(_spimage_pybackend ${SPIMAGE_LIBRARIES} "-undefined dynamic_lookup")
  else(APPLE)
    TARGET_LINK_LIBRARIES(_spimage_pybackend ${SPIMAGE_LIBRARIES} ${PYTHON_LIBRARIES})
  endif()

  set_target_properties(
    _spimage_pybackend
    PROPERTIES SOVERSION 1
    VERSION 1
    PREFIX ""
    INSTALL_NAME_DIR ${PYTHON_INSTDIR}
    )

  if(WIN32)
    set_target_properties(
      _spimage_pybackend
      PROPERTIES
      SUFFIX ".pyd"
      )
  else(WIN32)
    set_target_properties(
      _spimage_pybackend
      PROPERTIES
      SUFFIX ".so"
      )
  endif()
    
  INSTALL(TARGETS _spimage_pybackend
    RUNTIME DESTINATION ${PYTHON_INSTDIR}
    LIBRARY DESTINATION ${PYTHON_INSTDIR}
    ARCHIVE DESTINATION ${PYTHON_INSTDIR}
    )
  INSTALL(FILES spimage.py ${CMAKE_CURRENT_BINARY_DIR}/spimage_pybackend.py _spimage_conventions.py _spimage_cxi.py _spimage_resample.py _spimage_reconstructor.py _spimage_prtf.py _spimage_io.py _spimage_find_center.py _spimage_sphere_model.py _spimage_spheroid_model.py _spimage_pixelmask.py _spimage_qmap.py _spimage_utils.py _spimage_misc.py _spimage_patterson.py _spimage_material.py elements.dat leastsqbound-scipy/leastsqbound.py DESTINATION  ${PYTHON_INSTDIR})
else(SWIG_FOUND AND Python3_Interpreter_FOUND AND PYTHON_NUMPY_FOUND)
  MESSAGE("Skipping python wrappers: SWIG_FOUND=${SWIG_FOUND} Python3_Interpreter_FOUND=${Python3_Interpreter_FOUND} PYTHON_NUMPY_FOUND=${PYTHON_NUMPY_FOUND}")
endif()
