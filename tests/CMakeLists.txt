list(APPEND OBJECTS ${SPIMAGE_SRC})

if(CUDA_FOUND)
#  add_definitions(-D_STRICT_IEEE_754)
  list(APPEND OBJECTS "../src/fft.cu" "../src/phasing.cu" "../src/phasing_kernels.cu" "../src/image_filter_cuda.cu" "../src/support_update_cuda.cu")
endif()

if(WIN32)
set (GETOPT getopt.c)
add_definitions(-DGSL_DLL)
endif()

add_executable(tests ${OBJECTS} CuTest.c AllTests.c ImageTests.c LinearTests.c ProjTests.c ContainerTests.c PrtfTests.c PhasingTests.c ${GETOPT})
if(CUDA_FOUND)
  list(APPEND TESTS_LIBRARIES CUDA::curand CUDA::cufft)
  #set_target_properties(tests PROPERTIES CUDA_ARCHITECTURES "35")
endif()

target_link_libraries(tests PRIVATE ${TESTS_LIBRARIES})

enable_testing()

add_test(linearAlgebra tests -l)
add_test(phasing tests -p)
add_test(image tests -i)
add_test(projection tests -r)
add_test(PRTF tests -t)
add_test(container tests -c)
