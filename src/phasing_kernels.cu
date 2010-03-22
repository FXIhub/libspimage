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

__global__ void CUDA_support_projection_er(cufftComplex* g1,const int * pixel_flags,const  int size)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if((pixel_flags[i] & SpPixelInsideSupport) == 0){
      g1[i].x = 0;
      g1[i].y = 0;
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
  if(constraints & SpRampObject) {
    printf("Start ramp constraint\n");
    float x2 = 0.0;
    float y2 = 0.0;
    float xy = 0.0;
    float kx = 0.0;
    float ky = 0.0;
    float x_tmp;
    float y_tmp;
    float ax;
    float ay;
    for (int x = 0; x < sp_image_x(new_model); x++) {
      for (int y = 0; y < sp_image_y(new_model); y++) {
	if (sp_i3matrix_get(ph->pixel_flags,x,y,0) & SpPixelInsideSupport) {
	  if (x < sp_image_x(new_model)/2) {
	    x_tmp = (real) x;
	  } else {
	    x_tmp = (real)( sp_image_x(new_model) - x );
	  }
	  if (y < sp_image_y(new_model)/2) {
	    y_tmp = (real) y;
	  } else {
	    y_tmp = (real)( sp_image_y(new_model) - y );
	  }
	  x2 += x*x;
	  y2 += y*y;
	  xy += x*y;
	  //kx += x*sp_imag(sp_image_get(new_model,x,y,0));
	  //ky += y*sp_imag(sp_image_get(new_model,x,y,0));
	  kx += x*sp_carg(sp_image_get(new_model,x,y,0));
	  ky += y*sp_carg(sp_image_get(new_model,x,y,0));
	}
      }
    }
    ax = (kx*y2-ky*xy) / (x2*y2 - xy*xy);
    ay = (ky*x2-kx*xy) / (x2*y2 - xy*xy);
    for (int x = 0; x < sp_image_x(new_model); x++) {
      for (int y = 0; y < sp_image_y(new_model); y++) {
	if (sp_i3matrix_get(ph->pixel_flags,x,y,0) & SpPixelInsideSupport) {
	  if (x < sp_image_x(new_model)/2) {
	    x_tmp = (real) x;
	  } else {
	    x_tmp = (real)( sp_image_x(new_model) - x );
	  }
	  if (y < sp_image_y(new_model)/2) {
	    y_tmp = (real) y;
	  } else {
	    y_tmp = (real)( sp_image_y(new_model) - y );
	  }
	  sp_image_set(new_model,x,y,0,
		       sp_cinit(sp_cabs(sp_image_get(new_model,x,y,0))*
				cos(ax*x_tmp+ay*y_tmp),
				sp_cabs(sp_image_get(new_model,x,y,0))*
				sin(ax*x_tmp+ay*y_tmp)));
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
