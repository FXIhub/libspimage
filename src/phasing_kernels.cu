#include <cuda.h>
#include <cuda_runtime.h>
#include <spimage.h>

__global__ void CUDA_diff_map_f1(cufftComplex* f1, const cufftComplex* g0,const int * pixel_flags,const float gamma1,const  int size){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    f1[i] = g0[i];
    if((pixel_flags[i] & SpPixelInsideSupport) == 0){
      f1[i].x = -gamma1*g0[i].x;
      f1[i].y = -gamma1*g0[i].y;
    }
  }
}

__global__ void CUDA_diff_map(cufftComplex* Pi2f1,cufftComplex* Pi2rho, const cufftComplex* g0, cufftComplex* g1,const int * pixel_flags,const float gamma2,const float beta,const  int size){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if((pixel_flags[i] & SpPixelInsideSupport)){
      g1[i].x = g0[i].x +(beta)*((1+gamma2)*Pi2rho[i].x-gamma2*g0[i].x);
      g1[i].y = g0[i].y +(beta)*((1+gamma2)*Pi2rho[i].y-gamma2*g0[i].y);
    }else{
      g1[i] = g0[i];
    }
    g1[i].x -= beta*Pi2f1[i].x;
    g1[i].y -= beta*Pi2f1[i].y;
  }
}

__global__ void CUDA_support_projection_raar(cufftComplex* g1, const cufftComplex* g0, const cufftComplex* gp, const int * pixel_flags, const  int size, const float beta)
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
    if (pixel_flags[i] & SpPixelInsideSupport) {
      g1[i].x = gp[i].x;
      g1[i].y = gp[i].y;
    } else {
      g1[i].x = g0[i].x*beta+(1.0f-2.0f*beta)*gp[i].x;
      g1[i].y = g0[i].y*beta+(1.0f-2.0f*beta)*gp[i].y;
    }
  }
}      

__global__ void CUDA_support_projection_hio(cufftComplex* g1, const cufftComplex* g0, const cufftComplex* gp, const int * pixel_flags, const  int size, const float beta)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if (pixel_flags[i] & SpPixelInsideSupport) {
      g1[i].x = gp[i].x;
      g1[i].y = gp[i].y;
    } else {
      g1[i].x = g0[i].x-gp[i].x*beta;
      g1[i].y = g0[i].y-gp[i].y*beta;
    }
  }
}      

__global__ void CUDA_support_projection_er(cufftComplex* g1, cufftComplex *gp, const int * pixel_flags, const  int size)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(pixel_flags[i] & SpPixelInsideSupport) {
      g1[i].x = gp[i].x;
      g1[i].y = gp[i].y;
    } else {
      g1[i].x = 0;
      g1[i].y = 0;
    }
  }
}      


__global__ void CUDA_module_projection(cufftComplex* g, const float* amp, const float* amperrtol, const int * pixel_flags,const  int size, const SpPhasingConstraints constraints)
{	
  const int i =  blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(pixel_flags[i] & SpPixelMeasuredAmplitude){
      float m = 1.;
      if(!(constraints & SpAmplitudeErrorMargin) || (amperrtol==NULL)){
	// Default: Projection on measured amplitude
#ifndef _STRICT_IEEE_754      
        m = amp[i]/sqrt(g[i].x*g[i].x + g[i].y*g[i].y);     
#else
        m = __fdiv_rn(amp[i],__fsqrt_rn(__fadd_rn(__fmul_rn(g[i].x,g[i].x),
							    __fmul_rn( g[i].y,g[i].y))));     
#endif
      }else{
        // Projection according to given amplitude error tolerance map
#ifndef _STRICT_IEEE_754      
        const float amp_a = sqrt(g[i].x*g[i].x + g[i].y*g[i].y);     
#else
        const float amp_a = __fsqrt_rn(__fadd_rn(__fmul_rn(g[i].x,g[i].x),__fmul_rn( g[i].y,g[i].y)));     
#endif
        const float ampdiff = amp_a - amp[i];
        if (fabs(ampdiff) > amperrtol[i]){
          if (ampdiff > amperrtol[i]){
#ifndef _STRICT_IEEE_754      
            m = (amp[i]-amperrtol[i])/amp_a;
#else
            m = __fdiv_rn(amp[i]-amperrtol[i],amp_a);
#endif
	  }else{
#ifndef _STRICT_IEEE_754      
            m = (amp[i]+amperrtol[i])/amp_a;
#else
            m = __fdiv_rn(amp[i]+amperrtol[i],amp_a)
#endif
	  }
        }
      }
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

__global__ void CUDA_apply_fourier_constraints(cufftComplex* g,const  int size,const SpPhasingConstraints constraints){
  const int i =  blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(constraints & SpCentrosymmetricObject){
      if(g[i].x > 0){
	g[i].x  = sqrt(g[i].x*g[i].x+g[i].y*g[i].y);
      }else{
	g[i].x  = -sqrt(g[i].x*g[i].x+g[i].y*g[i].y);
      }
      g[i].y = 0;
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
