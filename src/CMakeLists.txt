LIST(APPEND OBJECTS  "image_util.c" "fft.c" "linear_alg.c" "mem_util.c" "image_sphere.c" "sperror.c" 
  "gaussianinv.c" "image_noise.c" "hashtable.c" "interpolation_kernels.c" "time_util.c" "list.c" "prtf.c")
IF(BUILD_STATIC_LIB)
	ADD_LIBRARY(spimage STATIC ${OBJECTS})
ELSE(BUILD_STATIC_LIB)
	ADD_LIBRARY(spimage SHARED ${OBJECTS})
ENDIF(BUILD_STATIC_LIB)


TARGET_LINK_LIBRARIES(spimage ${SPIMAGE_LIBRARIES})

INSTALL(TARGETS spimage
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIBRARY_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIBRARY_SUFFIX})
