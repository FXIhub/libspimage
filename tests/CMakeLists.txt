LIST(APPEND OBJECTS  "../src/image_util.c" "../src/fft.c" "../src/linear_alg.c" "../src/mem_util.c" "../src/image_sphere.c" "../src/sperror.c" 
  "../src/gaussianinv.c" "../src/image_noise.c" "../src/hashtable.c" "../src/interpolation_kernels.c" "../src/time_util.c" "../src/list.c" "../src/prtf.c" "../src/phasing.c" "../src/colormap.c" "../src/cuda_util.c" "../src/support_update.c")

IF(CUDA_FOUND)
#  ADD_DEFINITIONS(-D_STRICT_IEEE_754)
  LIST(APPEND OBJECTS "../src/fft.cu" "../src/phasing.cu" "../src/phasing_kernels.cu" "../src/image_filter_cuda.cu" "../src/support_update_cuda.cu")
ENDIF(CUDA_FOUND)


IF(CUDA_FOUND)
CUDA_ADD_EXECUTABLE(tests ${OBJECTS} CuTest.c AllTests.c ImageTests.c LinearTests.c ProjTests.c ContainerTests.c PrtfTests.c PhasingTests.c)
CUDA_ADD_CUFFT_TO_TARGET(tests)
ELSE(CUDA_FOUND)
ADD_EXECUTABLE(tests ${OBJECTS} CuTest.c AllTests.c ImageTests.c LinearTests.c ProjTests.c ContainerTests.c PrtfTests.c PhasingTests.c)
ENDIF(CUDA_FOUND)

TARGET_LINK_LIBRARIES(tests ${TESTS_LIBRARIES})

ENABLE_TESTING()

ADD_TEST(tests tests)