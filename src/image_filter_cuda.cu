#include "spimage.h"
//#include "cudpp/cudpp.h"


#include <thrust/sort.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <cstdlib>

static __global__ void CUDA_create_gaussian_kernel(cufftComplex * a,const int x,const int y, const int z,const float radius);
static __global__ void CUDA_Complex_multiply(cufftComplex * a,cufftComplex * b, int size);
static __global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale);

struct cmpCufftComplex{   
  __device__ bool operator()(const cufftComplex lhs, const cufftComplex rhs) { 
#ifdef _STRICT_IEEE_754
    if(__fadd_rn(__fmul_rn(lhs.x,lhs.x), __fmul_rn(lhs.y,lhs.y))
       > __fadd_rn(__fmul_rn(rhs.x,rhs.x), __fmul_rn(rhs.y,rhs.y))){
#else      
    if(lhs.x*lhs.x + lhs.y*lhs.y < rhs.x*rhs.x + rhs.y*rhs.y){
#endif
      return true;
    }
    return false;
  } 
};

struct addCufftComplex{   
  __device__ cufftComplex operator()(const cufftComplex lhs, const cufftComplex rhs) { 
    cufftComplex temp = lhs;
#ifdef _STRICT_IEEE_754
    temp.x = __fadd_rn(temp.x,rhs.x);
    temp.y = __fadd_rn(temp.y,rhs.y);
#else
    temp.x += rhs.x;
    temp.y += rhs.y;
#endif
    return temp;
  } 
};


void sp_gaussian_blur_cuda(cufftComplex * in, cufftComplex * out, int x, int y, int z, float radius, cufftHandle plan){
  cufftComplex * kernel = NULL;
  cutilSafeCall(cudaMalloc((void**)&kernel, sizeof(cufftComplex)*x*y*z));
  sp_create_gaussian_kernel_cuda(kernel,x,y,z,radius);
  cufftSafeCall(cufftExecC2C(plan, kernel, kernel, CUFFT_FORWARD));    
  cufftSafeCall(cufftExecC2C(plan, in, out, CUFFT_FORWARD));
  int blockSize = 0;
  int gridSize = 0;
  sp_cuda_launch_parameters(x*y*z, &gridSize, &blockSize);
  CUDA_Complex_multiply<<<gridSize,blockSize>>>(out,kernel,x*y*z);
  sp_cuda_check_errors();
  cufftSafeCall(cufftExecC2C(plan, out, out, CUFFT_INVERSE));
  CUDA_complex_scale<<<gridSize,blockSize>>>(out,x*y*z,1.0f/(x*y*z));
  sp_cuda_check_errors();
  cutilSafeCall(cudaFree(kernel));
}

// Complex pointwise multiplication
static __global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
#ifdef _STRICT_IEEE_754
    a[i].x = __fmul_rn(a[i].x,scale);
    a[i].y = __fmul_rn(a[i].y,scale);
#else
    a[i].x *= scale;
    a[i].y *= scale;
#endif
  }
} 

static __global__ void CUDA_Complex_multiply(cufftComplex * a,cufftComplex * b, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  float tmp_a_x;
  if(i<size){
    tmp_a_x = a[i].x;
#ifdef _STRICT_IEEE_754
    a[i].x = __fadd_rn(__fmul_rn(tmp_a_x,b[i].x),__fmul_rn(-a[i].y,b[i].y));
    a[i].y = __fadd_rn(__fmul_rn(tmp_a_x,b[i].y),__fmul_rn(a[i].y,b[i].x));
#else
    a[i].x = tmp_a_x*b[i].x-a[i].y*b[i].y;
    a[i].y = tmp_a_x*b[i].y+a[i].y*b[i].x;
#endif
  }
  
} 

static __global__ void CUDA_create_gaussian_kernel(cufftComplex * a,const int x,const int y, const int z,const float radius){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<x*y*z){
    const int my_x = i%x;
    const int my_y = (i/x)%y;
    const int my_z = i/(x*y);
    const int delta_z = min(my_z,z-my_z);
    const int delta_y = min(my_y,y-my_y);
    const int delta_x = min(my_x,x-my_x);
#ifdef _STRICT_IEEE_754
    const float d2 = __fadd_rn(__fadd_rn(__fmul_rn(delta_z,delta_z),__fmul_rn(delta_y,delta_y)),__fmul_rn(delta_x,delta_x));
    /* no IEEE strict exp unfortunately */
    const float f = exp(__fdiv_rn(-d2,(__fmul_rn(2,__fmul_rn(radius,radius)))));
#else
    const float d2 = delta_z*delta_z+delta_y*delta_y+delta_x*delta_x;
    const float f = exp(-d2/(2*radius*radius));
#endif
    a[my_z*(x*y)+my_y*x+my_x].x = f;
    a[my_z*(x*y)+my_y*x+my_x].y = 0;
  }
}

/*  The filter function is given by:

f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2)/(2*radius^2)) */
void sp_create_gaussian_kernel_cuda(cufftComplex * a,int x, int y, int z, float radius){
  int blockSize = 0;
  int gridSize = 0;
  sp_cuda_launch_parameters(x*y*z, &gridSize, &blockSize);
  // FM: This function is very slow for 3D images!!!
  CUDA_create_gaussian_kernel<<<gridSize,blockSize>>>(a,x,y,z,radius);
  sp_cuda_check_errors();
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(a);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(a+x*y*z));
  cufftComplex sum = {0,0};
  sum = thrust::reduce(beginc,endc,sum,addCufftComplex());
  sp_cuda_check_errors();
  sp_cuda_launch_parameters(x*y*z, &gridSize, &blockSize);
  CUDA_complex_scale<<<gridSize,blockSize>>>(a,x*y*z,1.0f/(sum.x+sum.y));
  sp_cuda_check_errors();
}

