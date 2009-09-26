#include "spimage.h"

static __global__ void CUDA_create_gaussian_kernel(cufftComplex * a,const int x,const int y, const int z,const float radius);
static __global__ void CUDA_Complex_multiply(cufftComplex * a,cufftComplex * b, int size);
static __global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale);
static __global__ void CUDA_slow_integrate(cufftComplex * a, int size ,float * result);

void sp_gaussian_blur_cuda(cufftComplex * in, cufftComplex * out, int x, int y, int z, float radius, cufftHandle plan){
  cufftComplex * kernel; 
  cutilSafeCall(cudaMalloc((void**)&kernel, sizeof(cufftComplex)*x*y*z));
  sp_create_gaussian_kernel_cuda(kernel,x,y,z,radius);
  cufftSafeCall(cufftExecC2C(plan, in, out, CUFFT_FORWARD));
  cufftSafeCall(cufftExecC2C(plan, kernel, kernel, CUFFT_FORWARD));
  int blockSize = 256;
  int gridSize = (x*y*z+blockSize-1)/blockSize;
  CUDA_Complex_multiply<<<gridSize,blockSize>>>(out,kernel,x*y*z);
  sp_cuda_check_errors();
  cufftSafeCall(cufftExecC2C(plan, out, out, CUFFT_INVERSE));
  CUDA_complex_scale<<<gridSize,blockSize>>>(out,x*y*z,1.0f/(x*y*z));
  sp_cuda_check_errors();
  cutilSafeCall(cudaFree(kernel));
  
}

// slow integrate
static __global__ void CUDA_slow_integrate(cufftComplex * a, int size ,float * result){
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  float reg = 0;
  if(id == 1){
    reg = 0;
    for(int i = 0;i<size;i++){
      reg += a[i].x;
    }
  }  
  *result = reg;
} 


// Complex pointwise multiplication
static __global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    a[i].x *= scale;
    a[i].y *= scale;
  }
} 

static __global__ void CUDA_Complex_multiply(cufftComplex * a,cufftComplex * b, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    a[i].x = a[i].x*b[i].x-a[i].y*b[i].y;
    a[i].y = a[i].x*b[i].y+a[i].y*b[i].x;
  }
  
} 

static __global__ void CUDA_create_gaussian_kernel(cufftComplex * a,const int x,const int y, const int z,const float radius){
  const int my_x = blockIdx.x*blockDim.x + threadIdx.x;
  const int my_y = blockIdx.y*blockDim.y + threadIdx.y;
  const int my_z = blockIdx.z*blockDim.z + threadIdx.z;
  if(my_x < x && my_y < y && my_z < z){
    const int delta_z = min(my_z,z-my_z);
    const int delta_y = min(my_y,y-my_y);
    const int delta_x = min(my_x,x-my_x);
    const float d2 = delta_z*delta_z+delta_y*delta_y+delta_x*delta_x;
    const float f = exp(-d2/(2*radius*radius));
    a[my_z*(x*y)+my_y*x+my_x].x = f;
    a[my_z*(x*y)+my_y*x+my_x].y = 0;
  }
}

/*  The filter function is given by:

f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2)/(2*radius^2)) */
void sp_create_gaussian_kernel_cuda(cufftComplex * a,int x, int y, int z, float radius){
  dim3 dimBlock(16,16,1);
  dim3 dimGrid((x+dimBlock.x-1)/dimBlock.x,
	       (y+dimBlock.y-1)/dimBlock.y,
	       (z+dimBlock.z-1)/dimBlock.z);
  CUDA_create_gaussian_kernel<<<dimGrid,dimBlock>>>(a,x,y,z,radius);
  sp_cuda_check_errors();
  float * d_sum;
  int blockSize = 256;
  int gridSize = (x*y*z+blockSize-1)/blockSize;
  cutilSafeCall(cudaMalloc((void**)&d_sum, sizeof(float)));
  CUDA_slow_integrate<<<gridSize,blockSize>>>( a, x*y*z ,d_sum);    
  sp_cuda_check_errors();
  float sum;
  cutilSafeCall(cudaMemcpy(&sum, d_sum,sizeof(float),cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaFree(d_sum));
  CUDA_complex_scale<<<gridSize,blockSize>>>(a,x*y*z,1.0f/sum);
  sp_cuda_check_errors();
}
