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

real sp_colormap_scale_value(real input,int colormap,real max_v, real min_v, double gamma){
  real scale,offset;
  if(min_v < 0){
    min_v = 0;
  }
  if(max_v < 0){
    max_v = 0;
  }
  if(max_v-min_v){
    scale = 65535/(max_v-min_v);
  }else{
    scale = 1;
  }
  offset = min_v;
  const real log_of_scale = log(65536);
  real value = sp_min(input,max_v);
  value = sp_max(value,min_v);
  value -= offset;
  value *= scale;
      
  if(colormap & SpColormapLogScale){
    value = log(value+1)/log_of_scale;
  }else{
    value /= 65535;
  }
  if(gamma != 1){
    value = pow(value,gamma);
  }
  return value;
}


void sp_colormap_write_rgb(unsigned char * out,const Image * img, int colormap,sp_rgb * color_table,real max_v, real min_v, int x, int y, int z,int red_blue_swap, double gamma){
  Complex cvalue = sp_image_get(img,x,y,z);
  real value = sp_colormap_scale_value(sp_cabs(cvalue),colormap,max_v,min_v,gamma);
  if(!isfinite(value)){   
    return;
  }
  if(colormap & SpColormapPhase){
    real phase = (256*(sp_carg(cvalue)+3.1416)/(2*3.1416));
    out[0] =  sqrt(value)*color_table[(int)phase].r;
    out[1] = sqrt(value)*color_table[(int)phase].g;
    out[2] = sqrt(value)*color_table[(int)phase].b;
  }else if(colormap & SpColormapMask){
    value = sp_image_mask_get(img,x,y,z);
    if(value){
      value = 255;
    }
    if(value){
      out[0] = color_table[(int)value].r;
      out[1] = color_table[(int)value].g;
      out[2] = color_table[(int)value].b;
    }else{
      /* use a checkered pattern to indicate no mask */
      if(((x%16)/8 + (y%16)/8) != 1){
	out[0] = 0x99;
	out[1] = 0x99;
	out[2] = 0x99;
      }else{
	out[0] = 0x66;
	out[1] = 0x66;
	out[2] = 0x66;
      }
    }    
  }else if(colormap & SpColormapShadedMask){    
    value *= 255;
    if(sp_image_mask_get(img,x,y,z)){
      out[0] = color_table[(int)value].r;
      out[1] = color_table[(int)value].g;
      out[2] = color_table[(int)value].b;
    }else{
      out[0] = color_table[(int)value].r/2;
      out[1] = color_table[(int)value].g/2;
      out[2] = color_table[(int)value].b/2;
    }
  }else{
    value *= 255;
    out[0] =  color_table[(int)value].r;
    out[1] = color_table[(int)value].g;
    out[2] = color_table[(int)value].b;    
  }  
  if(red_blue_swap){
    unsigned char tmp = out[0];
    out[0] = out[2];
    out[2] = tmp;
  }
}

sp_rgb sp_colormap_rgb_from_value(real value, int colormap){
  sp_rgb ret = {0,0,0};
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
  }else if(colormap & SpColormapBlackJet){
    if(value < 1/8.0){
      ret.r = 0;
      ret.g = 0;
      ret.b = (value)*8.0/3.0*4;
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
  }else if(colormap & SpColormapWhiteJet){
    if(value < 1/8.0){
      ret.r = (1/8.0-value)*8;
      ret.g = (1/8.0-value)*8;
      ret.b = 1;
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
  }else if(colormap & SpColormapWhiteJetShadow){
    if(value < 1/8.0){
      ret.r = (1/8.0-value)*8;
      ret.g = (1/8.0-value)*8;
      ret.b = (1/8.0-value)*8;
    }else if(value < 2/8.0){
      ret.r = 0;
      ret.g = 0;
      ret.b = (value-1.0/8.0)*8;
    }else if(value < 3/8.0){
      ret.r = 0;
      ret.g = (value-2.0/8.0)*8;
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
  }else if(colormap & SpColormapBlGrRdYl){
    float colors[][3] = {{0,0,0},{0,0,0},{0,0,35},{0,2,73},{4,64,103},{17,118,99},{25,151,58},{24,144,0},{98,108,0},{187,15,0},{206,43,0},{215,84,0},{225,125,0},{234,168,0},{244,211,0},{254,251,0}};
    float classes = sizeof(colors)/(3*sizeof(float));
    for(int i = 1;i<classes;i++){
      int j = i-1;
      float t = i/(classes-1);
      float u = (i-1)/(classes-1);
      if(value <= t){
	ret.r = colors[j][0]/255 + 2*(colors[i][0]-colors[j][0])*(value-u)/255.0;
	ret.g = colors[j][1]/255 + 2*(colors[i][1]-colors[j][1])*(value-u)/255.0;
	ret.b = colors[j][2]/255 + 2*(colors[i][2]-colors[j][2])*(value-u)/255.0;
	break;
      }
    }
  }else if(colormap & SpColormapCXI){
    float colors[][3] = {{0,0,0},{28,30,31},{0,24,255},{0,139,211},{247,134,0},{255,255,0},{255,255,255}};
    float classes = sizeof(colors)/(3*sizeof(float));
    for(int i = 1;i<classes;i++){
      int j = i-1;
      float t = i/(classes-1);
      float u = (i-1)/(classes-1);
      u = (value-u) * (classes-1);
      if(value <= t){
	ret.r = colors[j][0]/255 + (colors[i][0]-colors[j][0])*(u)/255.0;
	ret.g = colors[j][1]/255 + (colors[i][1]-colors[j][1])*(u)/255.0;
	ret.b = colors[j][2]/255 + (colors[i][2]-colors[j][2])*(u)/255.0;
	break;
      }
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

unsigned char * sp_image_get_false_color(Image * img, int color, double min, double max, double gamma){

  int x,y,z;
  sp_rgb color_table[256];
  real max_v,min_v;
  unsigned char * out = sp_malloc(sizeof(unsigned char)*sp_image_x(img)*sp_image_y(img)*sp_image_z(img)*4);

  max_v = max;
  min_v = min;

  sp_colormap_create_table(color_table,color);

  if(min == max){
    /* We're gonna scale the image so that it fits on the 8 bits */
    min_v = sp_c3matrix_min(img->image,NULL);
    max_v = sp_c3matrix_max(img->image,NULL);
  }
  for(z = 0;z<sp_image_z(img);z++){
    for(y = 0;y<sp_image_y(img);y++){
      for(x = 0;x<sp_image_x(img);x++){
	sp_colormap_write_rgb(&(out[4*(x+sp_image_x(img)*y+
				       z*sp_image_x(img)*
				       sp_image_y(img))]),img,color,
			      color_table,max_v,min_v,x,y,0,1,gamma);
	
      }
    }
  }
  return out;
}
