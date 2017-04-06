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
  cufftComplex * kernel; 
  cutilSafeCall(cudaMalloc((void**)&kernel, sizeof(cufftComplex)*x*y*z));
  sp_create_gaussian_kernel_cuda(kernel,x,y,z,radius);
  cufftSafeCall(cufftExecC2C(plan, in, out, CUFFT_FORWARD));
  cufftSafeCall(cufftExecC2C(plan, kernel, kernel, CUFFT_FORWARD));
  int blockSize = 512;
  int gridSize = (x*y*z+blockSize-1)/blockSize;
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
  const int my_x = blockIdx.x*blockDim.x + threadIdx.x;
  const int my_y = blockIdx.y*blockDim.y + threadIdx.y;
  /*  
      const int my_z = blockIdx.z*blockDim.z + threadIdx.z;
  */
  /* we're gonna have to do all z due to hardware limitations */
  for(int my_z = 0;my_z<z;my_z++){
    if(my_x < x && my_y < y && my_z < z){
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
}

/*  The filter function is given by:

f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2)/(2*radius^2)) */
void sp_create_gaussian_kernel_cuda(cufftComplex * a,int x, int y, int z, float radius){
  dim3 dimBlock(16,16,1);
  /* We have to use a grid size of 1 on z due to hardware limitations*/
  dim3 dimGrid((x+dimBlock.x-1)/dimBlock.x,
	       (y+dimBlock.y-1)/dimBlock.y,
	       1);
  CUDA_create_gaussian_kernel<<<dimGrid,dimBlock>>>(a,x,y,z,radius);
  sp_cuda_check_errors();
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(a);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(a+x*y*z));
  cufftComplex sum = {0,0};
  sum = thrust::reduce(beginc,endc,sum,addCufftComplex());
  sp_cuda_check_errors();
  int blockSize = 512;
  int gridSize = (x*y*z+blockSize-1)/blockSize;
  CUDA_complex_scale<<<gridSize,blockSize>>>(a,x*y*z,1.0f/(sum.x+sum.y));
  sp_cuda_check_errors();
}

/*
static __global__ void CUDA_complex_sum_reduce(cufftComplex * g_idata, cufftComplex * g_odata, unsigned int n){
  const int blockSize = blockDim.x;
  extern __shared__ cufftComplex sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockSize*2) + tid;
  unsigned int gridSize = blockSize*2*gridDim.x;
  sdata[tid].x = 0;
  sdata[tid].y = 0;
  
  while(i<n){
    sdata[tid].x += g_idata[i].x + g_idata[i+blockSize].x;
    sdata[tid].y += g_idata[i].y + g_idata[i+blockSize].y;
    i += gridSize;
  }
  __syncthreads();
  if(blockSize >= 512){
    if(tid < 256){
      sdata[tid].x += sdata[tid+256].x;
      sdata[tid].y += sdata[tid+256].y;
    }
    __syncthreads();
  }
  if(blockSize >= 256){
    if(tid < 128){
      sdata[tid].x += sdata[tid+128].x;
      sdata[tid].y += sdata[tid+128].y;
    }
    __syncthreads();
  }
  if(blockSize >= 128){
    if(tid < 64){
      sdata[tid].x += sdata[tid+64].x;
      sdata[tid].y += sdata[tid+64].y;
    }
    __syncthreads();
  }
  if(tid < 32){
    if(blockSize >= 64){
      sdata[tid].x += sdata[tid+32].x;
      sdata[tid].y += sdata[tid+64].y;
    }
    if(blockSize >= 32){
      sdata[tid].x += sdata[tid+16].x;
      sdata[tid].y += sdata[tid+16].y;
    }
    if(blockSize >= 16){
      sdata[tid].x += sdata[tid+8].x;
      sdata[tid].y += sdata[tid+8].y;
    }
    if(blockSize >= 8){
      sdata[tid].x += sdata[tid+4].x;
      sdata[tid].y += sdata[tid+4].y;
    }
    if(blockSize >= 4){
      sdata[tid].x += sdata[tid+2].x;
      sdata[tid].y += sdata[tid+2].y;
    }
    if(blockSize >= 2){
      sdata[tid].x += sdata[tid+1].x;
      sdata[tid].y += sdata[tid+1].y;
    }
  }
  if(tid == 0){
    g_odata[blockIdx.x].x = sdata[0].x;
    g_odata[blockIdx.x].y = sdata[0].y;
  }
}
*/
