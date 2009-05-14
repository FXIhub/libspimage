LIST(APPEND OBJECTS  "../src/image_util.c" "../src/fft.c" "../src/linear_alg.c" "../src/mem_util.c" "../src/image_sphere.c" "../src/sperror.c" 
  "../src/gaussianinv.c" "../src/image_noise.c" "../src/hashtable.c" "../src/interpolation_kernels.c" "../src/time_util.c" "../src/list.c" "../src/prtf.c" "../src/phasing.c")

ADD_EXECUTABLE(tests ${OBJECTS} CuTest.c AllTests.c ImageTests.c LinearTests.c ProjTests.c ContainerTests.c PrtfTests.c PhasingTests.c)

TARGET_LINK_LIBRARIES(tests ${TESTS_LIBRARIES})

ENABLE_TESTING()

ADD_TEST(tests tests)