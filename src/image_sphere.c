#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <tiffio.h>
#include <png.h>
#include <float.h>
#include <ctype.h>
#include <strings.h>
#include "spimage.h"

/*! This function calculates a rotation from three euler angles.
 */
SpRotation * sp_rot_alloc()
{
  SpRotation * rot = sp_matrix_alloc(3,3);
  return rot;
}

void sp_rot_free(SpRotation * rot)
{
  sp_matrix_free(rot);
  free(rot);
}

SpRotation * sp_rot_euler(real a1, real a2, real a3)
{
  sp_matrix * R1 = sp_matrix_alloc(3,3);
  sp_matrix * R2 = sp_matrix_alloc(3,3);
  sp_matrix * R3 = sp_matrix_alloc(3,3);
  sp_matrix_set(R1,0,0,cos(a1));
  sp_matrix_set(R1,1,1,cos(a1));
  sp_matrix_set(R1,0,1,sin(a1));
  sp_matrix_set(R1,1,0,-sin(a1));
  sp_matrix_set(R1,2,2,1.0);
  sp_matrix_set(R2,1,1,cos(a2));
  sp_matrix_set(R2,2,2,cos(a2));
  sp_matrix_set(R2,1,2,sin(a2));
  sp_matrix_set(R2,2,1,-sin(a2));
  sp_matrix_set(R2,0,0,1.0);
  sp_matrix_set(R3,0,0,cos(a3));
  sp_matrix_set(R3,1,1,cos(a3));
  sp_matrix_set(R3,0,1,sin(a3));
  sp_matrix_set(R3,1,0,-sin(a3));
  sp_matrix_set(R3,2,2,1.0);
  
  SpRotation * rot = sp_rot_alloc();
  sp_matrix_free(rot);

  sp_matrix * R4 = sp_matrix_mul(R2,R1);
  rot = sp_matrix_mul(R3,R4);

  sp_matrix_free(R1);
  sp_matrix_free(R2);
  sp_matrix_free(R3);
  sp_matrix_free(R4);
 
  return rot;
}


/* Taken from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToEuler/index.htm
 * might not work very well
 */
void sp_rot_get_euler(SpRotation * rot, real *alpha, real *beta, real *gamma)
{
  if (sp_matrix_get(rot,2,2) == 1.0) {
    *beta = M_PI;
    *alpha = atan2(sp_matrix_get(rot,0,1),sp_matrix_get(rot,1,1));
    *gamma = 0;
  } else if (sp_matrix_get(rot,2,2) == -1.0) {
    *beta = M_PI;
    *alpha = atan2(sp_matrix_get(rot,1,0),sp_matrix_get(rot,0,0));
    *gamma = 0;
  } else {
    *alpha = atan2(sp_matrix_get(rot,2,0),-sp_matrix_get(rot,2,1));
    *beta = acos(sp_matrix_get(rot,2,2));
    *gamma = atan2(sp_matrix_get(rot,0,2),sp_matrix_get(rot,1,2));
  }
}

SpRotation * sp_rot_multiply(SpRotation * a, SpRotation * b)
{
  SpRotation * rot = sp_rot_alloc();
  sp_matrix_free(rot);
  rot = sp_matrix_mul(a,b);
  return rot;
}

SpRotation * sp_rot_transpose(SpRotation * a)
{
  SpRotation * b = sp_rot_alloc();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      sp_matrix_set(b,i,j,sp_matrix_get(a,j,i));
  return b;
}

SpRotation * sp_rot_inverse()
{
  int i;
  SpRotation * rot = sp_rot_alloc();
  for (i = 0; i < 3; i++) {
    sp_matrix_set(rot,i,i,-1.0);
  }
  return rot;
}

/* Adds a disturbance to a rotation a. The disturbance is described as a small rotation along
 * the two degrees of freedom. They are both normally distributed with a standard deviation
 * of sigma (in radians)
 */
SpRotation * sp_rot_disturb(SpRotation * a, real sigma)
{
  int i;
  real theta = 0; //theta is normally idstributed
  real kappa = 0; //kappa is normally idstributed
  real phi = p_drand48()*2*M_PI;
  //this loop gives theta and kappa their (approximate) normal distribution
  for (i = 0; i < 12; i++){ theta += p_drand48(); kappa += p_drand48();}
  theta -= 6; theta *= sigma; kappa -= 6; kappa *= sigma;

  // R1 rotates the z-axis a little bit in the direction of phi.
  sp_matrix * R1 = sp_matrix_alloc(3,3);
  if (cos(theta) != 0) sp_matrix_set(R1,0,0,1.0/cos(theta));
  else sp_matrix_set(R1,2,0,1.0);
  sp_matrix_set(R1,1,1,1.0);
  sp_matrix_set(R1,0,2,cos(theta));
  sp_matrix_set(R1,1,2,sin(theta)*cos(phi));
  sp_matrix_set(R1,2,2,sin(theta)*sin(phi));

  // R2 rotates around the z-axis
  sp_matrix * R2 = sp_matrix_alloc(3,3);
  sp_matrix_set(R2,0,0,cos(kappa));
  sp_matrix_set(R2,1,1,cos(kappa));
  sp_matrix_set(R2,0,1,sin(kappa));
  sp_matrix_set(R2,1,0,-sin(kappa));

  SpRotation * rot = sp_rot_alloc();
  sp_matrix_free(rot);
  rot = sp_matrix_mul(sp_matrix_mul(R1,R2),a);
  return rot;
}

real sp_rot_difference(SpRotation * a, SpRotation * b)
{
  SpRotation * diff = sp_rot_multiply(a,sp_rot_transpose(b));
  real ret = (sp_matrix_get(diff,0,0)+sp_matrix_get(diff,1,1)+sp_matrix_get(diff,2,2))/3.0;
  sp_rot_free(diff);
  return ret;
}

real sp_rot_determinant(SpRotation * rot)
{
  return
    sp_matrix_get(rot,0,0)*sp_matrix_get(rot,1,1)*sp_matrix_get(rot,2,2)+
    sp_matrix_get(rot,1,0)*sp_matrix_get(rot,2,1)*sp_matrix_get(rot,0,2)+
    sp_matrix_get(rot,2,0)*sp_matrix_get(rot,0,1)*sp_matrix_get(rot,1,2)-
    sp_matrix_get(rot,0,0)*sp_matrix_get(rot,2,1)*sp_matrix_get(rot,1,2)-
    sp_matrix_get(rot,1,0)*sp_matrix_get(rot,0,1)*sp_matrix_get(rot,2,2)-
    sp_matrix_get(rot,2,0)*sp_matrix_get(rot,1,1)*sp_matrix_get(rot,0,2);
}

void sp_rot_draw(SpRotation * rot)
{
  printf("%f %f %f\n%f %f %f\n%f %f %f\n",
	 sp_matrix_get(rot,0,0),
	 sp_matrix_get(rot,0,1),
	 sp_matrix_get(rot,0,2),
	 sp_matrix_get(rot,1,0),
	 sp_matrix_get(rot,1,1),
	 sp_matrix_get(rot,1,2),
	 sp_matrix_get(rot,2,0),
	 sp_matrix_get(rot,2,1),
	 sp_matrix_get(rot,2,2));
}

SpRotation * sp_rot_uniform()
{
  real x1, x2, x3;
  sp_matrix * R = sp_matrix_alloc(3,3);
  sp_matrix * H = sp_matrix_alloc(3,3);
  x1 = p_drand48();
  x2 = p_drand48();
  x3 = p_drand48();

  sp_matrix_set(R,0,0,cos(2*M_PI*x1));
  sp_matrix_set(R,1,1,cos(2*M_PI*x1));
  sp_matrix_set(R,0,1,sin(2*M_PI*x1));
  sp_matrix_set(R,1,0,-sin(2*M_PI*x1));
  sp_matrix_set(R,2,2,1);

  sp_matrix_set(H,0,0,2.0*cos(2.0*M_PI*x2)*cos(2.0*M_PI*x2)*x3 - 1.0);
  sp_matrix_set(H,0,1,2.0*cos(2.0*M_PI*x2)*sin(2.0*M_PI*x2)*x3);
  sp_matrix_set(H,0,2,2.0*cos(2.0*M_PI*x2)*sqrt(x3-x3*x3));
  sp_matrix_set(H,1,0,2.0*cos(2.0*M_PI*x2)*sin(2.0*M_PI*x2)*x3);
  sp_matrix_set(H,1,1,2.0*sin(2.0*M_PI*x2)*sin(2.0*M_PI*x2)*x3 - 1.0);
  sp_matrix_set(H,1,2,2.0*sin(2.0*M_PI*x2)*sqrt(x3-x3*x3));
  sp_matrix_set(H,2,0,2.0*cos(2.0*M_PI*x2)*sqrt(x3-x3*x3));
  sp_matrix_set(H,2,1,2.0*sin(2.0*M_PI*x2)*sqrt(x3-x3*x3));
  sp_matrix_set(H,2,2,2.0*(1.0-x3) - 1.0);
  
   SpRotation * rot = sp_matrix_mul(H,R);
  
  sp_matrix_free(R);
  sp_matrix_free(H);

  return rot;
}

sp_3matrix * sp_image_sphere_z(Image * img)
{
  if (sp_image_z(img) != 1) {
    fprintf(stderr,"Trying to get z coordinates for a 3D image");
  }
  int x,y;
  real thisz;
  real D = img->detector->detector_distance;
  real px = img->detector->pixel_size[0];
  real py = img->detector->pixel_size[1];
  real cx = img->detector->image_center[0];
  real cy = img->detector->image_center[1];
  printf("z - alloc\n");
  sp_3matrix * z = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),1);
  printf("z - alloced\n");
  for(x=0; x<sp_3matrix_x(z); x++){
    for(y=0; y<sp_3matrix_y(z); y++){
      thisz = D-sqrt(D*D-((real)x-cx)*((real)x-cx)*px*px-((real)y-cy)*((real)y-cy)*py*py);
      sp_3matrix_set(z,x,y,0,thisz);
    }
  }
  printf("set values\n");
  if(z){
    return z;
  }else{
    fprintf(stderr,"Problem in creating z");
    exit(0);
  }
}

void sp_image_insert_ewald(Image * img, sp_3matrix * weight, Image * slice, sp_3matrix * curvature, SpRotation * rot, int kernel, real sigma, int radius) {
  int x,y;
  int xk,yk,zk;
  real nx,ny,nz;
  int nx_int, ny_int, nz_int;
  real w;
  real cx = slice->detector->image_center[0];
  real cy = slice->detector->image_center[1];
  real px = slice->detector->pixel_size[0];
  real py = slice->detector->pixel_size[1];

  for (x = 0; x < sp_3matrix_x(curvature); x++) {
    for (y = 0; y < sp_3matrix_y(curvature); y++) {
      //printf("-----------\n");
      //printf("%g, %g, %g\n",(real)x-cx,(real)y-cy,sp_3matrix_get(curvature,x,y,0));
      //printf("(%g %g %g)\n",sp_matrix_get(rot,0,0),sp_matrix_get(rot,1,0),sp_matrix_get(rot,2,0));
    
      //printf("%g %g %g\n",((real)x-cx)*px*sp_matrix_get(rot,0,0) + ((real)y-cy)*py*sp_matrix_get(rot,1,0) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,0),
      //((real)x-cx)*px*sp_matrix_get(rot,0,1) + ((real)y-cy)*py*sp_matrix_get(rot,1,1) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,1),
      //((real)x-cx)*px*sp_matrix_get(rot,0,2) + ((real)y-cy)*py*sp_matrix_get(rot,1,2) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,2));
      nx = (((real)x-cx)*px*sp_matrix_get(rot,0,0) + ((real)y-cy)*py*sp_matrix_get(rot,1,0) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,0))/
	img->detector->pixel_size[0] + img->detector->image_center[0];
      ny = (((real)x-cx)*px*sp_matrix_get(rot,0,1) + ((real)y-cy)*py*sp_matrix_get(rot,1,1) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,1))/
	img->detector->pixel_size[1] + img->detector->image_center[1];
      nz = (((real)x-cx)*px*sp_matrix_get(rot,0,2) + ((real)y-cy)*py*sp_matrix_get(rot,1,2) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,2))/
	img->detector->pixel_size[2] + img->detector->image_center[2];

      nx_int = (int)(nx+0.5);
      ny_int = (int)(ny+0.5);
      nz_int = (int)(nz+0.5);

      //printf("nx ny nz = %g %g %g\n",nx,ny,nz);
      if (kernel == SP_QSPLINE_KERNEL) radius = 2;
      if (kernel == SP_CSPLINE_KERNEL) radius = 3;
      if (kernel == SP_4SPLINE_KERNEL) radius = 3;

      for (xk = nx_int - radius; xk <= nx_int + radius; xk++) {
	for (yk = ny_int - radius; yk <= ny_int + radius; yk++) {
	  for (zk = nz_int - radius; zk <= nz_int + radius; zk++) {
	    if (xk >= 0 && xk < sp_image_x(img) && yk >= 0 && yk < sp_image_y(img) && zk >= 0 && zk < sp_image_z(img)) {
	      if (kernel == SP_NEAREST_KERNEL) {
		if (xk == nx_int && yk == ny_int && zk == nz_int) {
		  sp_image_set(img,xk,yk,zk,sp_cadd(sp_image_get(slice,x,y,0),sp_image_get(img,xk,yk,zk)));
		  sp_3matrix_set(weight,xk,yk,zk,sp_3matrix_get(weight,xk,yk,zk) + 1);
		}
	      }
	      if (kernel == SP_GAUSSIAN_KERNEL) {
		real r2 = ((real)xk - nx)*((real)xk - nx) + ((real)yk - ny)*((real)yk - ny) + ((real)zk - nz)*((real)zk - nz);
		if (r2 <= (real)(radius*radius)) {
		  w = exp(-r2/2/sigma/sigma)/sigma/sigma/sigma/15.7496; // (2*pi)^(3/2) = 15.7496
		  sp_image_set(img,xk,yk,zk,sp_cadd(sp_cscale(sp_image_get(slice,x,y,0),w),sp_image_get(img,xk,yk,zk))); 
		  sp_3matrix_set(weight,xk,yk,zk,sp_3matrix_get(weight,xk,yk,zk) + w);
		}
	      }
	      if (kernel == SP_LINEAR_KERNEL) {
		real r2 = ((real)xk - nx)*((real)xk - nx) + ((real)yk - ny)*((real)yk - ny) + ((real)zk - nz)*((real)zk - nz);
		if (r2 < (real)(radius*radius)) {
		  w = 1.0 - sqrt(r2)/(real)radius;
		  sp_image_set(img,xk,yk,zk,sp_cadd(sp_cscale(sp_image_get(slice,x,y,0),w),sp_image_get(img,xk,yk,zk))); 
		  sp_3matrix_set(weight,xk,yk,zk,sp_3matrix_get(weight,xk,yk,zk) + w);
		}
	      }
	      if (kernel == SP_QSPLINE_KERNEL) {
		real r = sqrt(((real)xk - nx)*((real)xk - nx) + ((real)yk - ny)*((real)yk - ny) + ((real)zk - nz)*((real)zk - nz));
		if (r < 1.5) {
		  if (r <= 0.5) w = -2.0*r*r + 1.5;
		  else w = r*r - 3.0*r + 2.25;
		  sp_image_set(img,xk,yk,zk,sp_cadd(sp_cscale(sp_image_get(slice,x,y,0),w),sp_image_get(img,xk,yk,zk))); 
		  sp_3matrix_set(weight,xk,yk,zk,sp_3matrix_get(weight,xk,yk,zk) + w);
		}
	      }
	      if (kernel == SP_CSPLINE_KERNEL) {
		real r = sqrt(((real)xk - nx)*((real)xk - nx) + ((real)yk - ny)*((real)yk - ny) + ((real)zk - nz)*((real)zk - nz));
		if (r < 2.0) {
		  if (r <= 1) w = 3.0*r*r*r - 2.0*r*r + 4.0;
		  else w = -1.0*r*r*r + 2.0*r*r + 4.0*r + 8.0;
		  sp_image_set(img,xk,yk,zk,sp_cadd(sp_cscale(sp_image_get(slice,x,y,0),w),sp_image_get(img,xk,yk,zk))); 
		  sp_3matrix_set(weight,xk,yk,zk,sp_3matrix_get(weight,xk,yk,zk) + w);
		}
	      }
	      if (kernel == SP_4SPLINE_KERNEL) {
		real r = sqrt(((real)xk - nx)*((real)xk - nx) + ((real)yk - ny)*((real)yk - ny) + ((real)zk - nz)*((real)zk - nz));
		if (r < 2.5) {
		  if (r <= 0.5) w = 6.0*r*r*r*r - 15.0*r*r + 14.375;
		  else if (r <= 1.5) w = -4.0*r*r*r*r + 20*r*r*r - 30.0*r*r + 5.0*r + 13.75;
		  else w = r*r*r*r - 10.0*r*r*r + 37.5*r*r - 62.5*r + 39.0625;
		  sp_image_set(img,xk,yk,zk,sp_cadd(sp_cscale(sp_image_get(slice,x,y,0),w),sp_image_get(img,xk,yk,zk)));
		  sp_3matrix_set(weight,xk,yk,zk,sp_3matrix_get(weight,xk,yk,zk) + w);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
}

Image * sp_image_generate_pattern(int side)
{
  int x,y,z;
  Image * real = sp_image_alloc(side,side,side);
  for(x=(int)(0.4*sp_image_x(real));x<(int)(0.6*sp_image_x(real));x++){
    for(y=(int)(0.4*sp_image_y(real));y<(int)(0.6*sp_image_y(real));y++){
      for(z=(int)(0.4*sp_image_z(real));z<(int)(0.6*sp_image_z(real));z++){
	if((x-side/2)*(x-side/2)+(y-side/2)*(y-side/2)/2+2*(z-side/2)*(z-side/2)<36.0){
	  sp_image_set(real,x,y,z,sp_cinit(10.5,0));
	}
      }
    }
  }
  sp_init_fft(1);
  Image * fourier = sp_image_fft(real);
  fourier->shifted = 1;
  sp_image_shift(fourier);
  return fourier;
}

void sp_image_get_2dpattern(Image * pattern, Image * slice, sp_3matrix * curvature,  SpRotation * rot, int kernel)
{
  int x,y;
  real nx,ny,nz;
  int nx_int,ny_int,nz_int;
  real cx = slice->detector->image_center[0];
  real cy = slice->detector->image_center[1];
  real px = slice->detector->pixel_size[0];
  real py = slice->detector->pixel_size[1];
  real w;
  real w_tot = 0;
  int x1,y1,z1;
  int radius;
  if (kernel == SP_QSPLINE_KERNEL){
    radius = 2;
  }else if (kernel == SP_CSPLINE_KERNEL){
    radius = 3;
  }else if (kernel == SP_4SPLINE_KERNEL){
    radius = 3;
  }else{
    sp_error_fatal("wrong kernel parameter given");
    return;
  }
  for (x = 0; x < sp_image_x(slice); x++) {
    for (y = 0; y < sp_image_y(slice); y++) {
      nx = (((real)x-cx)*px*sp_matrix_get(rot,0,0) + ((real)y-cy)*py*sp_matrix_get(rot,1,0) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,0))/
	pattern->detector->pixel_size[0] + pattern->detector->image_center[0];
      ny = (((real)x-cx)*px*sp_matrix_get(rot,0,1) + ((real)y-cy)*py*sp_matrix_get(rot,1,1) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,1))/
	pattern->detector->pixel_size[1] + pattern->detector->image_center[1];
      nz = (((real)x-cx)*px*sp_matrix_get(rot,0,2) + ((real)y-cy)*py*sp_matrix_get(rot,1,2) + sp_3matrix_get(curvature,x,y,0)*sp_matrix_get(rot,2,2))/
	pattern->detector->pixel_size[2] + pattern->detector->image_center[2];
      
      nx_int = (int)(nx+0.5);
      ny_int = (int)(ny+0.5);
      nz_int = (int)(nz+0.5);
      w_tot = 0;
      for (x1 = nx_int - radius; x1 <= nx_int + radius; x1++) {
	for (y1 = ny_int - radius; y1 <= ny_int + radius; y1++) {
	  for (z1 = nz_int - radius; z1 <= nz_int + radius; z1++) {
	    if (x1 >= 0 && x1 < sp_image_x(pattern) && y1 >= 0 && y1 < sp_image_y(pattern) && z1 >= 0 && z1 < sp_image_z(pattern)) {
	      if (kernel == SP_QSPLINE_KERNEL) {
		real r = sqrt(((real)x1 - nx)*((real)x1 - nx) + ((real)y1 - ny)*((real)y1 - ny) + ((real)z1 - nz)*((real)z1 - nz));
		if (r < 1.5) {
		  if (r <= 0.5) w = -2*r*r + 1.5;
		  else w = r*r - 3*r + 2.25;
		  sp_image_set(slice,x,y,0,sp_cadd(sp_cscale(sp_image_get(pattern,x1,y1,z1),w),sp_image_get(pattern,x1,y1,z1))); 
		  w_tot += w;
		}
	      }
	      if (kernel == SP_CSPLINE_KERNEL) {
		real r = sqrt(((real)x1 - nx)*((real)x1 - nx) + ((real)y1 - ny)*((real)y1 - ny) + ((real)z1 - nz)*((real)z1 - nz));
		if (r < 2.0) {
		  if (r <= 1.0) w = 3*r*r*r - 2*r*r + 4.0;
		  else w = -1*r*r*r + 2*r*r + 4*r + 8.0;
		  sp_image_set(slice,x,y,0,sp_cadd(sp_cscale(sp_image_get(pattern,x1,y1,z1),w),sp_image_get(pattern,x1,y1,z1))); 
		  w_tot += w;
		}
	      }
	      if (kernel == SP_4SPLINE_KERNEL) {
		real r = sqrt(((real)x1 - nx)*((real)x1 - nx) + ((real)y1 - ny)*((real)y1 - ny) + ((real)z1 - nz)*((real)z1 - nz));
		if (r < 2.5) {
		  if (r <= 0.5) w = 6.0*r*r*r*r - 15.0*r*r + 14.375;
		  else if (r <= 1.5) w = -4.0*r*r*r*r + 20*r*r*r - 30.0*r*r + 5.0*r + 13.75;
		  else w = r*r*r*r - 10.0*r*r*r + 37.5*r*r - 62.5*r + 39.0625;
		  sp_image_set(slice,x,y,0,sp_cadd(sp_cscale(sp_image_get(pattern,x1,y1,z1),w),sp_image_get(pattern,x1,y1,z1))); 
		  w_tot += w;
		}
	      }
	    }
	  }
	}
      }
    }
    if (w_tot > 0.01)
      sp_image_set(slice,x,y,0,sp_cscale(sp_image_get(slice,x,y,0),1/w_tot));
  }
}

// I want to make it possible to up/downsample when slicing. a is the downsampling factor
void sp_image_get_slice(Image * space, Image * slice, sp_3matrix * slice_z, real a, SpRotation * rot)
{
  int x,y;
  real nx,ny,nz;
  sp_c3matrix * local;
  
  real cx = slice->detector->image_center[0];
  real cy = slice->detector->image_center[1];
  real space_cx = space->detector->image_center[0];
  real space_cy = space->detector->image_center[1];
  real space_cz = space->detector->image_center[2];

  local = sp_c3matrix_alloc(sp_image_x(slice),sp_image_y(slice),sp_image_z(slice));

  for(x=0; x<sp_c3matrix_x(local); x++){
    for(y=0; y<sp_c3matrix_y(local); y++){
      nx = a*(x-cx)*sp_matrix_get(rot,0,0)+a*(y-cy)*sp_matrix_get(rot,1,0)+a*sp_3matrix_get(slice_z,x,y,0)/slice->detector->pixel_size[0]*sp_matrix_get(rot,2,0);
      ny = a*(x-cx)*sp_matrix_get(rot,0,1)+a*(y-cy)*sp_matrix_get(rot,1,1)+a*sp_3matrix_get(slice_z,x,y,0)/slice->detector->pixel_size[0]*sp_matrix_get(rot,2,1);
      nz = a*(x-cx)*sp_matrix_get(rot,0,2)+a*(y-cy)*sp_matrix_get(rot,1,2)+a*sp_3matrix_get(slice_z,x,y,0)/slice->detector->pixel_size[0]*sp_matrix_get(rot,2,2);
      if(nx+space_cx>=0 && nx+space_cx<sp_image_x(space) &&
	 ny+space_cy>=0 && ny+space_cy<sp_image_y(space) &&
	 nz+space_cz>=0 && nz+space_cz<sp_image_z(space)){
	/*  I believe that sp_image_interp is the best here but for some reason the 
	 *  result looks really bad with it.
	 */
	sp_image_set(slice,x,y,0,sp_c3matrix_get(space->image,(int)nx+space_cx,(int)ny+space_cy,(int)nz+space_cz));
	//sp_image_set(slice,x,y,0,sp_c3matrix_get(space->image,(int)(nx+space_cx+0.5),(int)(ny+space_cy+0.5),(int)(nz+space_cz)));
	//sp_image_set(slice,x,y,0,sp_c3matrix_interp(space->image,nx+space_cx,ny+space_cy,nz+space_cz));
	/*
	get_value(space,local,weight,x,y,(int)nx,(int)ny,(int)nz,(nx-(int)nx)*(ny-(int)ny)*(nz-(int)nz));
	get_value(space,local,weight,x,y,(int)nx,(int)ny,(int)nz+1,(nx-(int)nx)*(ny-(int)ny)*(1-(int)nz+nz));
	get_value(space,local,weight,x,y,(int)nx,(int)ny+1,(int)nz,(nx-(int)nx)*(1-(int)ny+ny)*(nz-(int)nz));
	get_value(space,local,weight,x,y,(int)nx+1,(int)ny,(int)nz,(1-(int)nx+nx)*(ny-(int)ny)*(nz-(int)nz));
	get_value(space,local,weight,x,y,(int)nx,(int)ny+1,(int)nz+1,(nx-(int)nx)*(1-(int)ny+ny)*(1-(int)nz+nz));
	get_value(space,local,weight,x,y,(int)nx+1,(int)ny,(int)nz+1,(1-(int)nx+nx)*(ny-(int)ny)*(1-(int)nz+nz));
	get_value(space,local,weight,x,y,(int)nx+1,(int)ny+1,(int)nz,(1-(int)nx+1)*(1-(int)ny+1)*(nz-(int)nz));
	get_value(space,local,weight,x,y,(int)nx+1,(int)ny+1,(int)nz+1,(1-(int)nx+nx)*(1-(int)ny+ny)*(1-(int)nz+nz));
	*/
      }else{
	sp_i3matrix_set(slice->mask,x,y,0,1);
      }
    }
  }
}

