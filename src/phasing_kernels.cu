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


__global__ void CUDA_module_projection(cufftComplex* g, const float* amp, const float* amp_min, const float* amp_max, const int * pixel_flags,const  int size, const SpPhasingConstraints constraints)
{	
  const int i =  blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(pixel_flags[i] & SpPixelMeasuredAmplitude){
      float m = 1.;
      if(!(constraints & SpAmplitudeErrorMargin) || (amp_min==NULL) || (amp_max==NULL)){
	// Default: Projection on measured amplitude
#ifndef _STRICT_IEEE_754      
        m = amp[i]/sqrt(g[i].x*g[i].x + g[i].y*g[i].y);     
#else
        m = __fdiv_rn(amp[i],__fsqrt_rn(__fadd_rn(__fmul_rn(g[i].x,g[i].x),
							    __fmul_rn( g[i].y,g[i].y))));     
#endif
      }else{
        // Projection according to given amplitude margins
#ifndef _STRICT_IEEE_754      
        const float amp_a = sqrt(g[i].x*g[i].x + g[i].y*g[i].y);     
#else
        const float amp_a = __fsqrt_rn(__fadd_rn(__fmul_rn(g[i].x,g[i].x),__fmul_rn( g[i].y,g[i].y)));     
#endif
	if (amp_a < amp_min[i]){
#ifndef _STRICT_IEEE_754      
          m = amp_min[i]/amp_a;
#else
          m = __fdiv_rn(amp_min,amp_a);
#endif
	}
	else if(amp_a > amp_max[i]){
#ifndef _STRICT_IEEE_754      
          m = amp_max[i]/amp_a;
#else
          m = __fdiv_rn(amp_max[i],amp_a)
#endif
	}
      }
      if(isfinite(m)){
#ifndef _STRICT_IEEE_754      
	g[i].x *= m;
	g[i].y *= m;
#else
	g[i].x = __fmul_rn(g[i].x,m);
	g[i].y = __fmul_rn(g[i].y,m);
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

// Complex pointwise absolute value squared
__global__ void CUDA_complex_abs2(cufftComplex * a, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    a[i].x = a[i].x*a[i].x+a[i].y*a[i].y;
    a[i].y = 0;
  }
}

// Complex pointwise real space error
__global__ void CUDA_ereal(cufftComplex * out, const cufftComplex * in, const int * pixel_flags, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    /* We'll store the denominator in the imaginary part */
    float abs2 = in[i].x*in[i].x+in[i].y*in[i].y;
    out[i].y = abs2;
    /* And the numerator in the real part */
    if(pixel_flags[i] & SpPixelInsideSupport){
      /* inside the support */
      out[i].x = 0;
    }else{
      out[i].x = abs2;
    }
  }
}

// Stores the support in the real part and the intensities mask in the imaginary part
__global__ void CUDA_pixel_flags_to_complex(cufftComplex * out, const int * pixel_flags, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    if(pixel_flags[i] & SpPixelInsideSupport){
      out[i].x = 1;
    }else{
      out[i].x = 0;
    }
    if(pixel_flags[i] & SpPixelMeasuredAmplitude){
      out[i].y = 1;
    }else{
      out[i].y = 0;
    }    
  }
}

// Complex pointwise fourier space error
__global__ void CUDA_efourier(cufftComplex * out, const cufftComplex * fmodel, const float* amp, const int * pixel_flags, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    /* We'll store the denominator in the imaginary part */
    /* And the numerator in the real part */    
    float fmodel_abs2 = fmodel[i].x*fmodel[i].x+fmodel[i].y*fmodel[i].y;
    if(pixel_flags[i] & SpPixelMeasuredAmplitude){
      out[i].y = amp[i]*amp[i];
      out[i].x = amp[i]-sqrt(fmodel_abs2);
      /* squared */
      out[i].x *= out[i].x;
    }else{
      out[i].y = fmodel_abs2;
      out[i].x = 0;
    }
  }
}

// Complex pointwise amplitude ratio
__global__ void CUDA_FcFo(cufftComplex * out, const cufftComplex * fmodel, const float* amp, const int * pixel_flags, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    /* We'll store the denominator in the imaginary part */
    /* And the numerator in the real part */    
    float fmodel_abs2 = fmodel[i].x*fmodel[i].x+fmodel[i].y*fmodel[i].y;
    if(pixel_flags[i] & SpPixelMeasuredAmplitude){
      out[i].y = amp[i];
      out[i].x = sqrt(fmodel_abs2);
    }else{
      out[i].y = 0;
      out[i].x = 0;
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

// Complex pointwise addition
__global__ void CUDA_complex_add(cufftComplex * a, int size ,cufftComplex add){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<size){
    a[i].x += add.x;
    a[i].y += add.y;
  }
}

// Add random phases while preserving amplitude
// The input uniform_random must have uniformly distributed number from 0 to 1
__global__ void CUDA_random_rephase(cufftComplex * a, float * uniform_random, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i<size){
    const float amp_a = sqrt(a[i].x*a[i].x + a[i].y*a[i].y);
    float phase = uniform_random[i]*2*M_PI;
    a[i].x = cos(phase)*amp_a;
    a[i].y = sin(phase)*amp_a;
  }
}


__global__ void CUDA_real_to_complex(cufftComplex * out, float * in, int size){
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i<size){
    out[i].x = in[i];
    out[i].y = 0;
  }
}
