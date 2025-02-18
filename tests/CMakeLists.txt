LIST(APPEND OBJECTS ${SPIMAGE_SRC})

IF(CUDA_FOUND)
#  ADD_DEFINITIONS(-D_STRICT_IEEE_754)
  LIST(APPEND OBJECTS "../src/fft.cu" "../src/phasing.cu" "../src/phasing_kernels.cu" "../src/image_filter_cuda.cu" "../src/support_update_cuda.cu")
ENDIF(CUDA_FOUND)

IF(WIN32)
SET (GETOPT getopt.c)
ADD_DEFINITIONS(-DGSL_DLL)
ENDIF(WIN32)

IF(CUDA_FOUND)
CUDA_ADD_EXECUTABLE(tests ${OBJECTS} CuTest.c AllTests.c ImageTests.c LinearTests.c ProjTests.c ContainerTests.c PrtfTests.c PhasingTests.c ${GETOPT})
CUDA_ADD_CUFFT_TO_TARGET(tests)
target_link_libraries(tests ${CUDA_curand_LIBRARY})
ELSE(CUDA_FOUND)
ADD_EXECUTABLE(tests ${OBJECTS} CuTest.c AllTests.c ImageTests.c LinearTests.c ProjTests.c ContainerTests.c PrtfTests.c PhasingTests.c ${GETOPT})
ENDIF(CUDA_FOUND)

TARGET_LINK_LIBRARIES(tests ${TESTS_LIBRARIES})

ENABLE_TESTING()

ADD_TEST(linearAlgebra tests -l)
ADD_TEST(phasing tests -p)
ADD_TEST(image tests -i)
ADD_TEST(projection tests -r)
ADD_TEST(PRTF tests -t)
ADD_TEST(container tests -c)
