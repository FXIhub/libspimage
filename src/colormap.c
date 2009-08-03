#include "spimage.h"

static void hsv_to_rgb(float H,float S,float V,float * R,float *G,float *B);


static void hsv_to_rgb(float H,float S,float V,float * R,float *G,float *B){
  if( V == 0 ){ 
    *R = 0;
    *G = 0;
    *B = 0; 
  }else if( S == 0 ) {                                                                   
    *R = V;                                                            
    *G = V;                                                            
    *B = V;                                                            
  } else {                                                                   
    const double hf = H / 60.0;                                       
    const int    i  = (int) floor( hf );                              
    const double f  = hf - i;                                         
    const double pv  = V * ( 1 - S );                                 
    const double qv  = V * ( 1 - S * f );                             
    const double tv  = V * ( 1 - S * ( 1 - f ) );                     
    switch( i ){                                                               
    case 0: 
      *R = V; 
      *G = tv;
      *B = pv;
      break; 
    case 1:
      *R = qv;
      *G = V;
      *B = pv;
      break;
    case 2:
      *R = pv; 
      *G = V;
      *B = tv;
      break; 
    case 3: 
      *R = pv;
      *G = qv;
      *B = V;
      break;
    case 4:  
      *R = tv; 
      *G = pv;
      *B = V;
      break;  
    case 5:
      *R = V;
      *G = pv;
      *B = qv; 
      break;
    case 6: 
      *R = V;
      *G = tv;    
      *B = pv; 
      break; 
    case -1:  
      *R = V;
      *G = pv; 
      *B = qv;
      break;
    default:
      sp_error_fatal("i Value error in HSV to *R*G*B conversion, Value is %d",i);
      break;
    }									
  }									
  *R *= 255.0F;                                                        
  *G *= 255.0F;                                                        
  *B *= 255.0F;  
}


sp_rgb sp_colormap_rgb_from_value(real value, int colormap){
  sp_rgb ret;
  if(colormap & SpColormapGrayScale){       
    ret.r = value;
    ret.g = value;
    ret.b = value;
  }else if(colormap & SpColormapTraditional){       
    ret.r = sqrt(value);
    ret.g = value*value*value;
    ret.b = sin(value*2*M_PI);
  }else if(colormap & SpColormapHot){
    ret.r = 3*value;
    ret.g = 3*value-1;
    ret.b = 3*value-2;
  }else if(colormap & SpColormapRainbow){
    ret.r = fabs(2*value-0.5);
    ret.g = sin(value*M_PI);
    ret.b = cos(value*M_PI/2);
  }else if(colormap & SpColormapJet){
    if(value < 1/8.0){
      ret.r = 0;
      ret.g = 0;
      ret.b = (value+1.0/8.0)*4;	   
    }else if(value < 3/8.0){
      ret.r = 0;
      ret.g = (value-1.0/8.0)*4;
      ret.b = 1;
    }else if(value < 5/8.0){
      ret.r = (value-3.0/8.0)*4;
      ret.g = 1;
      ret.b = 1-(value-3.0/8.0)*4;
    }else if(value < 7/8.0){
      ret.r = 1;
      ret.g = 1-(value-5.0/8.0)*4;
      ret.b = 0;
    }else if(value <= 1.01){
      ret.r = 1-(value-7.0/8.0)*4;;
      ret.g = 0;
      ret.b = 0;
    }
  }

  ret.r = sp_min(1,ret.r);
  ret.g = sp_min(1,ret.g);
  ret.b = sp_min(1,ret.b);
  ret.r = sp_max(0,ret.r);
  ret.g = sp_max(0,ret.g);
  ret.b = sp_max(0,ret.b);
  ret.r *= 255;
  ret.g *= 255;
  ret.b *= 255;

  if(colormap & SpColormapWheel){
    hsv_to_rgb(360*value,1.0,1.0,&ret.r,&ret.g,&ret.b);
  }
  return ret;
}


void sp_colormap_create_table(sp_rgb color_table[256],int colormap){
  for(int i = 0;i<256;i++){
    real value = i/255.0;
    color_table[i] = sp_colormap_rgb_from_value(value,colormap);
  }
}
