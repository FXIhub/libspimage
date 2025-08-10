list(APPEND OBJECTS ${SPIMAGE_SRC})



if(CUDA_FOUND)
  list(APPEND OBJECTS "fft.cu" "phasing.cu" "phasing_kernels.cu" "support_update_cuda.cu" "image_filter_cuda.cu")
endif()

if(BUILD_STATIC_LIB)
  add_library(spimage STATIC ${OBJECTS})
else()
add_library(spimage SHARED ${OBJECTS})
if(CUDA_FOUND)
  list(APPEND SPIMAGE_LIBRARIES CUDA::curand CUDA::cufft)
  #set_target_properties(spimage PROPERTIES CUDA_ARCHITECTURES "35")
endif()
endif()

target_link_libraries(spimage PRIVATE ${SPIMAGE_LIBRARIES})

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
endif()


install(TARGETS spimage
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(SWIG_FOUND AND Python3_Interpreter_FOUND AND PYTHON_NUMPY_FOUND)
  file(GLOB spimage_headers "${CMAKE_CURRENT_SOURCE_DIR}/../include/spimage/*.h")

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c COMMAND ${SWIG_EXECUTABLE} -python -outcurrentdir  -o ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c ${CMAKE_CURRENT_SOURCE_DIR}/spimage_pybackend.i DEPENDS ${spimage_headers} ${CMAKE_CURRENT_SOURCE_DIR}/spimage_pybackend.i)
  include_directories(${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR})

  add_library(_spimage_pybackend SHARED ${OBJECTS} ${CMAKE_CURRENT_BINARY_DIR}/spimage_wrap.c)

  if(APPLE)
    target_link_libraries(_spimage_pybackend ${SPIMAGE_LIBRARIES} "-undefined dynamic_lookup")
  else()
    target_link_libraries(_spimage_pybackend ${SPIMAGE_LIBRARIES} ${PYTHON_LIBRARIES})
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
  else()
    set_target_properties(
      _spimage_pybackend
      PROPERTIES
      SUFFIX ".so"
      )
  endif()
    
  install(TARGETS _spimage_pybackend
    RUNTIME DESTINATION ${PYTHON_INSTDIR}
    LIBRARY DESTINATION ${PYTHON_INSTDIR}
    ARCHIVE DESTINATION ${PYTHON_INSTDIR}
    )
  install(FILES spimage.py ${CMAKE_CURRENT_BINARY_DIR}/spimage_pybackend.py _spimage_conventions.py _spimage_cxi.py _spimage_resample.py _spimage_reconstructor.py _spimage_prtf.py _spimage_io.py _spimage_find_center.py _spimage_sphere_model.py _spimage_spheroid_model.py _spimage_pixelmask.py _spimage_qmap.py _spimage_utils.py _spimage_misc.py _spimage_patterson.py _spimage_material.py elements.dat leastsqbound-scipy/leastsqbound.py DESTINATION  ${PYTHON_INSTDIR})
else()
  message("Skipping python wrappers: SWIG_FOUND=${SWIG_FOUND} Python3_Interpreter_FOUND=${Python3_Interpreter_FOUND} PYTHON_NUMPY_FOUND=${PYTHON_NUMPY_FOUND}")
endif()
