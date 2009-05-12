#include "spimage.h"




Image * sp_prtf_basic(Image ** list,int n){
  if(n < 1){
    return NULL;
  }
  for(int i = 0;i<n;i++){
    if(sp_image_x(list[0]) != sp_image_x(list[i]) ||
       sp_image_y(list[0]) != sp_image_y(list[i]) ||
       sp_image_z(list[0]) != sp_image_z(list[i])){
      return NULL;
    }
       
  }
  Image * ret = sp_image_duplicate(list[0],SP_COPY_ALL);
  int size = sp_image_size(list[0]);
  sp_image_fill(ret,sp_cinit(0,0));
  for(int i = 0;i<n;i++){
    sp_image_add(ret,list[i]);
  }
  for(int i = 0;i<size;i++){
    real div = 0;
    for(int j = 0;j<n;j++){
      div += sp_cabs(list[j]->image->data[i]);
    }
    if(div){
      ret->image->data[i] = sp_cscale(ret->image->data[i],1.0/div);
    }
  }
  return ret;
}



Image * sp_prtf_advanced(Image ** list,int n,int flags){
  Image ** real_image = sp_malloc(sizeof(Image *)*n);
  Image ** fourier_image = sp_malloc(sizeof(Image *)*n);
  /* We'll do the superposition and phase match in real space
     but the prtf in reciprocal space
  */
  if((flags & (SpInPlace | SpOutOfPlace)) == 0){
    return NULL;
  }
  if(flags & SpRealSpace){
    for(int i = 0;i<n;i++){
      real_image[i] = sp_image_duplicate(list[i],SP_COPY_ALL);
    }
  }else if(flags & SpFourierSpace){
    for(int i = 0;i<n;i++){
      real_image[i] = sp_image_ifft(list[i]);
    }
  }else{
    return NULL;
  }
  for(int i = 1;i<n;i++){
    sp_image_superimpose_fractional(real_image[0],real_image[i],SpCorrectPhaseShift|SpEnantiomorph,4);
  }  
  for(int i = 0;i<n;i++){
    fourier_image[i] = sp_image_fft(real_image[i]);
    sp_image_scale(fourier_image[i],1.0/sp_image_size(fourier_image[i]));
  }
  Image * ret = sp_prtf_basic(fourier_image,n);
  for(int i = 0;i<n;i++){
    if(flags & SpInPlace){
      if(flags & SpFourierSpace){
	sp_image_memcpy(list[i],fourier_image[i]);
      }else{
	sp_image_memcpy(list[i],real_image[i]);
      }
    }
    sp_image_free(real_image[i]);
    sp_image_free(fourier_image[i]);
  }
  sp_free(fourier_image);
  sp_free(real_image);
  return ret;
}

sp_list * sp_prtf_by_resolution(Image * prtf){
  int size = sp_image_x(prtf)+sp_image_y(prtf)+sp_image_z(prtf);
  real * sum = sp_malloc(sizeof(real)*size);
  int * div = sp_malloc(sizeof(int)*size);
  for(int i = 0;i<size;i++){
    sum[i] = 0;
    div[i] = 0;
  }
  for(int i = 0;i<sp_image_size(prtf);i++){
    int d = sp_image_dist(prtf,i,SP_TO_CENTER);
    sum[d] += sp_cabs(prtf->image->data[i]);
    div[d]++;
  }
  int max_size = 0;
  for(int i = 0;i<size;i++){
    if(div[i]){
      max_size = i;
    }
  }
  sp_list * ret = sp_list_alloc(max_size);
  for(int i = 0;i<max_size;i++){
    if(div[i]){
      sp_list_append(ret,sum[i]/div[i]);
    }else{
      sp_list_append(ret,0);
    }
  }
  sp_free(sum);
  sp_free(div);
  return ret;
}
