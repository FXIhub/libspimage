#include <stdlib.h>
#include <string.h>
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif
#include "spimage.h"

sp_kernel * sp_create_spline2_kernel_table_by_size(real max_r2,int size){
  sp_kernel * ret = sp_malloc(sizeof(sp_kernel));
  ret->table = sp_malloc(sizeof(real)*size);
  ret->boundary = max_r2;
  ret->boundary_in_3D = (-3 + 2 * sqrt(3)*sqrt(max_r2))/6;
  ret->type = Kernel_Real_R2;
  ret->size = size;
  for(int i = 0;i<ret->size;i++){
    real r = sqrt(max_r2*i/ret->size);
    real w = 0;
    if(r < 1.5){
      if(r < 0.5){	    
	w = (r+1.5)*(r+1.5) - 3.0*(r+0.5)*(r+0.5);
      }else if(r < 1.5){
	w = (r+1.5)*(r+1.5) - 3.0*(r+0.5)*(r+0.5)+3.0*(r-0.5)*(r-0.5);
      }
    }
    ret->table[i] = w;
  }
  return ret;
}


sp_kernel * sp_create_spline2_kernel_table(real max_r2,real tolerance){
  /* According to Mathematica and after much thinking I think the error
   on the approximation is given by:
   -((3 (1/2 - Sqrt[1/4 + step] + step (5/2 + Sqrt[1/4 + step] - 2 Sqrt[step + 2 (1/4 + 1/2 Sqrt[1/4 + step])])))/(4 step))
   
   And the step size is given by:
   step -> 1/18*(36*tolerance + 16*tolerance*tolerance + sqrt(432*tolerance + pow(-36*tolerance - 16*tolerance*tolerance,2)))   
  */
  real step = 1.0/18*(36*tolerance + 16*tolerance*tolerance + sqrt(432*tolerance + pow(-36*tolerance - 16*tolerance*tolerance,2)));
  int size = (int)max_r2/step+1;
  return sp_create_spline2_kernel_table_by_size(max_r2,size);
}


sp_kernel * sp_create_spline2_r_kernel_table(real max_r,int size){
  sp_kernel * ret = sp_malloc(sizeof(sp_kernel));
  ret->table = sp_malloc(sizeof(real)*size);
  ret->boundary = max_r;
  ret->type = Kernel_Real_R;
  ret->size = size;
  for(int i = 0;i<ret->size;i++){
    real r = max_r*i/ret->size;
    real w = 0;
    if(r < 1.5){
      if(r < 0.5){	    
	w = (r+1.5)*(r+1.5) - 3.0*(r+0.5)*(r+0.5);
      }else if(r < 1.5){
	w = (r+1.5)*(r+1.5) - 3.0*(r+0.5)*(r+0.5)+3.0*(r-0.5)*(r-0.5);
      }
    }
    ret->table[i] = w;
  }
  return ret;
}

real sp_kernel_table_sample(sp_kernel * k, real r2){
  /* nearest neighbour interpolation */
  /*      return k->table[(int)(k->size*r2/k->boundary+0.5)]; */
  /* linear interpolation */
  real pos = k->size*r2/k->boundary;
  real u = pos-(int)pos;
  real ret = (1-u)*k->table[(int)(pos)]+(u)*k->table[(int)(pos)+1];
  return ret;
}

void sp_kernel_free(sp_kernel * k){
  sp_free(k->table);
  free(k);
}
