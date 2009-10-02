#include <cuda.h>
#include <cuda_runtime.h>
#include <spimage.h>

__global__ void CUDA_support_projection_raar(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta)
{
  /* A bit of documentation about the equation:
     
     Rs = 2*Ps-I; Rm = 2*Pm-I
     
     RAAR = 1/2 * beta * (RsRm + I) + (1 - beta) * Pm;    
     RAAR = 2*beta*Ps*Pm+(1-2*beta)*Pm - beta * (Ps-I)
     
     Which reduces to:
     
     Inside the support: Pm
     Outside the support: (1 - 2*beta)*Pm + beta*I
     
  */    
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if((pixel_flags[i] & SpPixelInsideSupport) == 0){
      g1[i].x = g0[i].x*beta+(1.0f-2.0f*beta)*g1[i].x;
      g1[i].y = g0[i].y*beta+(1.0f-2.0f*beta)*g1[i].y;
    }
  }
}      

__global__ void CUDA_support_projection_hio(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if((pixel_flags[i] & SpPixelInsideSupport) == 0){
      g1[i].x = g0[i].x-g1[i].x*beta;
      g1[i].y = g0[i].y-g1[i].y*beta;
    }
  }
}      

__global__ void CUDA_module_projection(cufftComplex* g, const float* amp,const int * pixel_flags,const  int size)
{	
  const int i =  blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(pixel_flags[i] & SpPixelMeasuredAmplitude){
#ifndef _STRICT_IEEE_754      
      const float m = amp[i]/sqrt(g[i].x*g[i].x + g[i].y*g[i].y);     
#else
      const float m = __fdiv_rn(amp[i],__fsqrt_rn(__fadd_rn(__fmul_rn(g[i].x,g[i].x),
							    __fmul_rn( g[i].y,g[i].y))));     
#endif
      if(isfinite(m)){
#ifndef _STRICT_IEEE_754      
	g[i].x *= m;
	g[i].y *= m;
#else
	g[i].x = __fmul_rn(g[i].x,m);
	g[i].y = __fmul_rn(g[i].y,m);;
#endif
      }else{
	g[i].x = amp[i];
	g[i].y = 0;
      }
    }
  }
}  


__global__ void CUDA_phased_amplitudes_projection(cufftComplex* g, const cufftComplex* phased_amplitudes,const int * pixel_flags,const  int size)
{	
  const int i =  blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(pixel_flags[i] & SpPixelMeasuredAmplitude){
      g[i].x = phased_amplitudes[i].x;
      g[i].y = phased_amplitudes[i].y;
    }
  }
}  


__global__ void CUDA_apply_constraints(cufftComplex* g, const int * pixel_flags,const  int size,const SpPhasingConstraints constraints){
  const int i =  blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(pixel_flags[i] & SpPixelInsideSupport){
      if(constraints & SpRealObject){
	g[i].y = 0;
      }else if(constraints & SpPositiveRealObject){
	if(g[i].x < 0){
	  if(constraints & SpPositivityFlipping){
	    g[i].x = fabs(g[i].x);
	  }else{
	    g[i].x = 0;
	  }
	}
	g[i].y = 0;
      }else if(constraints & SpPositiveComplexObject){
	if(g[i].x < 0){
	  if(constraints & SpPositivityFlipping){
	    g[i].x = fabs(g[i].x);
	  }else{
	    g[i].x = 0;
	  }
	}
	if(g[i].y < 0){
	  if(constraints & SpPositivityFlipping){
	    g[i].y = fabs(g[i].y);
	  }else{
	    g[i].y = 0;
	  }
	}
      }
    }
  }
}

// Complex pointwise multiplication
__global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    a[i].x *= scale;
    a[i].y *= scale;
  }
} 
