#ifndef IMAGETESTS_H_
#define IMAGETESTS_H_ 1

#include "AllTests.h"


void test_sp_image_shift(CuTest * tc){
  /* First test a shift on an even image without need for padding */
  Image * a = sp_image_alloc(2,2,1);
  a->detector->image_center[0] = 1;
  a->detector->image_center[1] = 1;
  a->detector->image_center[2] = 0;
  sp_image_set(a,0,0,0,sp_cinit(0,0));
  sp_image_set(a,1,0,0,sp_cinit(1,0));
  sp_image_set(a,0,1,0,sp_cinit(2,0));
  sp_image_set(a,1,1,0,sp_cinit(3,0));
  Image * b = sp_image_shift(a);
  CuAssertIntEquals(tc,sp_image_size(a),sp_image_size(b));
  CuAssertIntEquals(tc,sp_image_x(a),sp_image_x(b));
  CuAssertIntEquals(tc,sp_image_y(a),sp_image_y(b));
  CuAssertIntEquals(tc,sp_image_z(a),sp_image_z(b));
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,1,0,0),sp_image_get(a,0,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,1,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,0,1,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,1,1,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,1,1,0)))+REAL_EPSILON));  


  /* Now test to shift back */
  Image * c = sp_image_shift(b);
  CuAssertIntEquals(tc,sp_image_size(a),sp_image_size(c));
  CuAssertIntEquals(tc,sp_image_x(a),sp_image_x(c));
  CuAssertIntEquals(tc,sp_image_y(a),sp_image_y(c));
  CuAssertIntEquals(tc,sp_image_z(a),sp_image_z(c));
  CuAssertComplexEquals(tc,sp_image_get(c,0,0,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,0,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,1,0,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,1,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,0,1,0),sp_image_get(a,0,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,0,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,1,1,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,1,1,0)))+REAL_EPSILON));  

  sp_image_free(a);
  sp_image_free(b);
  sp_image_free(c);

  /* Now test an even image with need for padding*/
  a = sp_image_alloc(2,2,1);
  a->detector->image_center[0] = 0;
  a->detector->image_center[1] = 0;
  a->detector->image_center[2] = 0;
  sp_image_set(a,0,0,0,sp_cinit(0,0));
  sp_image_set(a,1,0,0,sp_cinit(1,0));
  sp_image_set(a,0,1,0,sp_cinit(2,0));
  sp_image_set(a,1,1,0,sp_cinit(3,0));

  b = sp_image_shift(a);
  /* This time we must have two zero paddings on the right size of the image*/
  CuAssertIntEquals(tc,4,sp_image_x(b));
  CuAssertIntEquals(tc,4,sp_image_y(b));
  CuAssertIntEquals(tc,1,sp_image_z(b));

  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,1,0,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,1,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,0,1,0),sp_image_get(a,0,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,1,1,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,1,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,2,1,0),sp_cinit(0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,2,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,1,2,0),sp_cinit(0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,2,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,2,2,0),sp_cinit(0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,2,1,0)))+REAL_EPSILON));  

  /* Now test to shift back */
  c = sp_image_shift(b);
  CuAssertIntEquals(tc,sp_image_size(b),sp_image_size(c));
  CuAssertIntEquals(tc,sp_image_x(b),sp_image_x(c));
  CuAssertIntEquals(tc,sp_image_y(b),sp_image_y(c));
  CuAssertIntEquals(tc,sp_image_z(b),sp_image_z(c));
  CuAssertComplexEquals(tc,sp_image_get(c,2,2,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,2,2,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,3,2,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,3,2,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,2,3,0),sp_image_get(a,0,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,2,3,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,3,3,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,3,3,0)))+REAL_EPSILON));    

  sp_image_free(a);
  sp_image_free(b);
  sp_image_free(c);

  /* Now test an odd image, with need for padding */

  a = sp_image_alloc(3,3,1);
  a->detector->image_center[0] = 1;
  a->detector->image_center[1] = 1;
  a->detector->image_center[2] = 0;
  sp_image_set(a,0,0,0,sp_cinit(0,0));
  sp_image_set(a,1,0,0,sp_cinit(1,0));
  sp_image_set(a,0,1,0,sp_cinit(2,0));
  sp_image_set(a,1,1,0,sp_cinit(3,0));

  b = sp_image_shift(a);
  /* This time we must have two zero paddings on the right size of the image*/
  CuAssertIntEquals(tc,4,sp_image_x(b));
  CuAssertIntEquals(tc,4,sp_image_y(b));
  CuAssertIntEquals(tc,1,sp_image_z(b));

  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,1,0,0),sp_image_get(a,2,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,1,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,0,1,0),sp_image_get(a,1,2,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,1,1,0),sp_image_get(a,2,2,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,1,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,3,3,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,1,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,3,0,0),sp_image_get(a,0,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,3,0,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(b,0,3,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,3,0)))+REAL_EPSILON));  


  /* Now test to shift back */
  c = sp_image_shift(b);
  CuAssertIntEquals(tc,sp_image_size(b),sp_image_size(c));
  CuAssertIntEquals(tc,sp_image_x(b),sp_image_x(c));
  CuAssertIntEquals(tc,sp_image_y(b),sp_image_y(c));
  CuAssertIntEquals(tc,sp_image_z(b),sp_image_z(c));
  CuAssertComplexEquals(tc,sp_image_get(c,1,1,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,1,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,2,1,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,2,1,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,1,2,0),sp_image_get(a,0,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,1,2,0)))+REAL_EPSILON));  
  CuAssertComplexEquals(tc,sp_image_get(c,2,2,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(c,2,2,0)))+REAL_EPSILON));    
  sp_image_free(a);
  sp_image_free(b);
  sp_image_free(c);

}

void test_image_print(Image * a){
  for(int z = 0; z<sp_image_z(a);z++){
    printf("***  Layer z = %d ***\n",z);
    for(int y = sp_image_y(a)-1; y>=0;y--){
      for(int x = 0; x<sp_image_x(a);x++){
	printf("(%3.2g,%3.2g)\t",sp_real(sp_image_get(a,x,y,z)),sp_imag(sp_image_get(a,x,y,z)));
      }
      printf("\n");
    }
    printf("\n");
  }
}

void test_sp_image_edge_extend(CuTest * tc){
  Image * a = sp_image_alloc(2,2,1);
  sp_image_set(a,0,0,0,sp_cinit(1,0));
  sp_image_set(a,1,0,0,sp_cinit(2,0));
  sp_image_set(a,0,1,0,sp_cinit(3,0));
  sp_image_set(a,1,1,0,sp_cinit(4,0));
  Image * b = sp_image_edge_extend(a,1,SP_ZERO_PAD_EDGE,SP_2D);
  for(int x = 0;x<sp_image_x(b);x++){
    for(int y = 0;y<sp_image_y(b);y++){
      if(!x || x == 3 || !y || y == 3){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,0),sp_cinit(0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,x,y,0)))+REAL_EPSILON));
      }
    }
  }
  sp_image_free(b);
  b = sp_image_edge_extend(a,1,SP_SYMMETRIC_EDGE,SP_2D);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,0,0)))+REAL_EPSILON));
  CuAssertComplexEquals(tc,sp_image_get(b,2,0,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,2,0,0))+REAL_EPSILON));
  sp_image_free(b);
  b = sp_image_edge_extend(a,1,SP_REPLICATE_EDGE,SP_2D);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,0,0,0))+REAL_EPSILON));
  CuAssertComplexEquals(tc,sp_image_get(b,2,0,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,2,0,0))+REAL_EPSILON));
  sp_image_free(b);
  b = sp_image_edge_extend(a,1,SP_CIRCULAR_EDGE,SP_2D);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,0,0,0))+REAL_EPSILON));
  CuAssertComplexEquals(tc,sp_image_get(b,2,0,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,2,0,0))+REAL_EPSILON));
  sp_image_free(b);
  sp_image_free(a);
}

void test_sp_bubble_sort(CuTest * tc){
  real array[4] = {4,-5.3, 1000,2};
  sp_bubble_sort(array,4);
  CuAssertDblEquals(tc,array[0],-5.3,fabs(REAL_EPSILON*(-5.3)+REAL_EPSILON));
  CuAssertDblEquals(tc,array[1],2,fabs(REAL_EPSILON*(2)+REAL_EPSILON));
  CuAssertDblEquals(tc,array[2],4,fabs(REAL_EPSILON*(4)+REAL_EPSILON));
  CuAssertDblEquals(tc,array[3],1000,fabs(REAL_EPSILON*(1000)+REAL_EPSILON));  
}


void test_sp_image_median_filter(CuTest * tc){
  Image * a = sp_image_alloc(2,2,1);
  sp_image_set(a,0,0,0,sp_cinit(1,0));
  sp_image_set(a,1,0,0,sp_cinit(2,0));
  sp_image_set(a,0,1,0,sp_cinit(3,0));
  sp_image_set(a,1,1,0,sp_cinit(4,0));
  sp_i3matrix * kernel = sp_i3matrix_alloc(2,2,1);
  sp_i3matrix_add_constant(kernel,2);
  sp_image_median_filter(a,kernel,SP_ZERO_PAD_EDGE,SP_2D);
  CuAssertDblEquals(tc,sp_cabs(sp_image_get(a,0,0,0)),0.0,fabs(REAL_EPSILON*(0)+REAL_EPSILON));
  CuAssertDblEquals(tc,sp_cabs(sp_image_get(a,1,1,0)),2.5,fabs(REAL_EPSILON*(2.5)+REAL_EPSILON));
  sp_image_free(a);
  sp_i3matrix_free(kernel);
}

void test_sp_image_gaussian_blur(CuTest * tc){
  Image * a = sp_image_alloc(11,11,1);
  sp_image_fill(a,sp_cinit(0,0));
  sp_image_set(a,5,5,0,sp_cinit(1,0));
  Image * b = gaussian_blur(a,1);
  /* Check if integral over image is mantained */
  Complex a_sum = {0,0};
  Complex b_sum = {0,0};
  for(int i = 0;i<sp_image_size(a);i++){
    sp_cincr(a_sum,a->image->data[i]);
    sp_cincr(b_sum,b->image->data[i]);
  }
  /* the 0.00015 follows from the fact that 99.7% of the gaussian is within 3
     standard deviation and the gaussian blur only calculates 3 standard deviations. */
  CuAssertComplexEquals(tc,a_sum,b_sum,fabs(0.00015*(sp_cabs(a_sum))));
  long long ind1;
  long long ind2;
  /* Image maximum must remain in the same place */
  sp_image_max(a,&ind1,NULL,NULL,NULL);
  sp_image_max(b,&ind2,NULL,NULL,NULL);
  CuAssertIntEquals(tc,ind1,ind2);
  sp_image_free(a);
  sp_image_free(b);
  a = sp_image_alloc(11,11,11);
  sp_image_fill(a,sp_cinit(0,0));
  /* sp_image_set(a,5,5,5,sp_cinit(1,0));
  sp_image_set(a,3,5,5,sp_cinit(1,0));
  sp_image_set(a,5,1,5,sp_cinit(1,0)); */
  sp_image_set(a,1,5,3,sp_cinit(1,0));
  b = gaussian_blur(a,1);
  sp_image_write(a,"3D_gaussian_blur_test_original.vtk",0);
  sp_image_write(b,"3D_gaussian_blur_test_blurred.vtk",0);
  /* Check if integral over image is mantained */
  a_sum = sp_cinit(0,0);
  b_sum = sp_cinit(0,0);
  for(int i = 0;i<sp_image_size(a);i++){
    sp_cincr(a_sum,a->image->data[i]);
    sp_cincr(b_sum,b->image->data[i]);
  }
  /* the 0.00015 follows from the fact that 99.7% of the gaussian is within 3
     standard deviation and the gaussian blur only calculates 3 standard deviations. */
  CuAssertComplexEquals(tc,a_sum,b_sum,fabs(0.00015*(sp_cabs(a_sum))));
  /* Image maximum must remain in the same place */
  sp_image_max(a,&ind1,NULL,NULL,NULL);
  sp_image_max(b,&ind2,NULL,NULL,NULL);
  CuAssertIntEquals(tc,ind1,ind2);
  sp_image_free(a);
  sp_image_free(b);
  
}

void test_cube_crop(CuTest * tc){
  Image * a = sp_image_alloc(3,3,3);
  for(long long i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(2,0);
  }
  sp_image_set(a,2,1,1,sp_cinit(3,0));
  Image *b = cube_crop(a,2,1,1,2,1,1);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_cinit(3,0),3*fabs(REAL_EPSILON));  
}


void test_sp_image_low_pass(CuTest * tc){
  Image * a = sp_image_alloc(3,3,4);
  Image * b;
  sp_image_set(a,0,0,0,sp_cinit(3.4,0));
  a->shifted = 1;
  b = sp_image_low_pass(a,1,0);  
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,0,0,0),3*fabs(REAL_EPSILON));  
  sp_image_free(b);
  sp_image_set(a,0,0,3,sp_cinit(3.1,0));
  b = sp_image_low_pass(a,2,0);  
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,2),sp_image_get(a,0,0,3),3*fabs(REAL_EPSILON));  
  
  sp_image_free(a);
  sp_image_free(b);
}
  


void test_sp_image_h5_read_write(CuTest * tc){
  /* First lets check 2D */
  int size = 100;
  Image * a;
  Image * b;
  a = sp_image_alloc(size,size,1);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }
  sp_image_write(a,"test.h5",sizeof(float));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");

  a = sp_image_alloc(size,size,1);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }

  /* Now with double precision */
  sp_image_write(a,"test.h5",sizeof(double));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");

  /* And now 3D */
  size = 10;
  a = sp_image_alloc(size,size,size);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }
  sp_image_write(a,"test.h5",sizeof(float));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");


  /* And in double precision */
  size = 10;
  a = sp_image_alloc(size,size,size);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }
  sp_image_write(a,"test.h5",sizeof(double));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");
}


void test_sp_image_get_false_color(CuTest * tc){
  int size = 100;
  Image * a;
  a = sp_image_alloc(size,size,1);
  a->phased = 1;
  unsigned char * out;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,x,y,z,sp_cinit(rand()%50,rand()%50));
      }
    }
  }
  out = sp_image_get_false_color(a,COLOR_JET,0,0);
  free(out);
  sp_image_get_false_color(a,COLOR_JET|LOG_SCALE,0,0);
  free(out);
  sp_image_get_false_color(a,COLOR_JET|LOG_SCALE,0,25);
  free(out);
}

void test_sp_image_noise_estimate(CuTest * tc){
  Image * object = sp_image_alloc(100,100,1);
  object->detector->image_center[0] = 50;
  object->detector->image_center[1] = 50;
  /* Create an circle with radius 10 */
  for(int i = 0;i<sp_image_size(object);i++){
    if(sp_image_dist(object,i,SP_TO_CENTER) <= 10){
      object->image->data[i] = sp_cinit(1,0);
    }else{
      object->image->data[i] = sp_cinit(0,0);
    }
  }
  sp_image_write(object,"test_noise_object.vtk",0);
  Image * p_intensities = sp_image_fft(object);
  p_intensities->scaled = 1;
  sp_image_to_intensities(p_intensities);
  sp_image_write(p_intensities,"test_noise_p_intensities.vtk",0);
  Image * p_autocorrelation = sp_image_ifft(p_intensities);
  sp_image_write(p_autocorrelation,"test_noise_p_autocorrelation.vtk",0);
  Image * autocorrelation_support = sp_image_alloc(sp_image_x(object),sp_image_y(object),sp_image_z(object));
  for(int i = 0;i<sp_image_size(object);i++){
    if(sp_real(p_autocorrelation->image->data[i]) < 1){
      autocorrelation_support->image->data[i] = sp_cinit(0,0);
    }else{
      autocorrelation_support->image->data[i] = sp_cinit(1,0);
    }
  }
  sp_image_write(autocorrelation_support,"test_noise_autocorrelation_support.vtk",0);

  Image * n_intensities = sp_image_alloc(sp_image_x(object),sp_image_y(object),sp_image_z(object));
  Image * p_std_dev = sp_image_alloc(sp_image_x(object),sp_image_y(object),sp_image_z(object));
  for(int i = 0;i<sp_image_size(object);i++){
    p_std_dev->image->data[i] = sp_cinit(sqrt(sp_real(p_intensities->image->data[i])/10),0);
    /* The noisy image will have points with standard deviation equal to the sqrt of the intensity */
    sp_real(n_intensities->image->data[i]) = sp_box_muller(sp_real(p_intensities->image->data[i]),sqrt(sp_real(p_intensities->image->data[i]))/10);
    sp_imag(n_intensities->image->data[i]) = 0;
    if(sp_real(n_intensities->image->data[i]) < 0){
      sp_real(n_intensities->image->data[i]) = 0;
    }
  }
  sp_image_write(n_intensities,"test_noise_n_intensities.vtk",0);
  sp_image_write(n_intensities,"test_noise_p_std_dev.vtk",0);
  Image * n_std_dev = sp_image_noise_estimate(n_intensities,autocorrelation_support);
  sp_image_write(n_std_dev,"test_noise_n_std_dev.vtk",0);
  
  Image * n_p_std_dev_ratio = sp_image_alloc(sp_image_x(object),sp_image_y(object),sp_image_z(object));
  for(int i = 0;i<sp_image_size(object);i++){
    sp_real(n_p_std_dev_ratio->image->data[i]) = sp_real(n_std_dev->image->data[i])/(sp_real(p_std_dev->image->data[i])+1);
    sp_imag(n_p_std_dev_ratio->image->data[i]) = 0;
  }
  sp_image_write(n_p_std_dev_ratio,"test_noise_n_p_std_dev_ratio.vtk",0);
  Image * noise = sp_image_alloc(sp_image_x(object),sp_image_y(object),sp_image_z(object));
  for(int i = 0;i<sp_image_size(object);i++){
    sp_real(noise->image->data[i]) = fabs(sp_real(n_intensities->image->data[i])-(sp_real(p_intensities->image->data[i])));
    sp_imag(noise->image->data[i]) = 0;
  }
  sp_image_write(noise,"test_noise_noise.vtk",0);
}

void test_sp_image_dist(CuTest * tc){
  int size = 100;
  Image * a;
  a = sp_image_alloc(size,size,1);
  /* Images default with the center in the center*/
  int i = sp_image_get_index(a,5,5,0);
  real dist = sp_image_dist(a,i,SP_TO_CORNER);
  CuAssertDblEquals(tc,dist,sqrt(2*5*5),100*dist*REAL_EPSILON);
   dist = sp_image_dist(a,i,SP_TO_CENTER);
  CuAssertDblEquals(tc,dist,sqrt(2*44.5*44.5),100*dist*REAL_EPSILON);
   dist = sp_image_dist(a,i,SP_TO_CENTER2);
  CuAssertDblEquals(tc,dist,44.5,100*dist*REAL_EPSILON);
   dist = sp_image_dist(a,i,SP_TO_AXIS);
  CuAssertDblEquals(tc,dist,0,100*dist*REAL_EPSILON);
}

void test_sp_image_max(CuTest * tc){
  int size = 100;
  Image * a;
  int x,y,z;
  long long index;
  a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(0,0);
  }
  sp_image_set(a,10,4,0,sp_cinit(2,4));
  sp_image_max(a,&index,&x,&y,&z);
  CuAssertIntEquals(tc,x,10);
  CuAssertIntEquals(tc,y,4);
  CuAssertIntEquals(tc,z,0);
  CuAssertIntEquals(tc,index,sp_image_get_index(a,x,y,z));
  sp_image_free(a);
}

void test_sp_image_superimpose(CuTest * tc){
  /* First a simple test without enantiomorph */
  int size = 3;
  Image * a;
  Image * b;
  a = sp_image_alloc(size,size,1);
  b = sp_image_alloc(size,size,1);
  a->phased = 1;
  b->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,x,y,z,sp_cinit(0,0));
	sp_image_set(b,x,y,z,sp_cinit(0,0));
      }
    }
  }  
  /*
         A       B
 
       0 1 0   0 0 0
       0 2 0   0 0 1
       0 0 0   0 0 2
   */
  sp_image_set(a,1,2,0,sp_cinit(1,1));
  sp_image_set(a,1,1,0,sp_cinit(2,-1));
  sp_image_set(b,1+1,2-1,0,sp_cinit(1,1));
  sp_image_set(b,1+1,1-1,0,sp_cinit(2,-1));
  sp_image_superimpose(a,b,0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }
  /* Now a test with enantiomorph */
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,x,y,z,sp_cinit(0,0));
	sp_image_set(b,x,y,z,sp_cinit(0,0));
      }
    }
  }  
  /*
         A       B
 
       0 1 0   0 0 0
       0 2 3   0 3 2
       0 0 0   0 0 1
   */
  sp_image_set(a,1,2,0,sp_cinit(1,1));
  sp_image_set(a,1,1,0,sp_cinit(2,-1));
  sp_image_set(a,2,1,0,sp_cinit(3,-1));

  sp_image_set(b,2,0,0,sp_cinit(1,1));
  sp_image_set(b,2,1,0,sp_cinit(2,-1));
  sp_image_set(b,1,1,0,sp_cinit(3,-1));

  sp_image_superimpose(a,b,SP_ENANTIOMORPH);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }
}

void test_sp_image_reflect(CuTest * tc){
  /* First a simple test without enantiomorph */
  int size = 3;
  Image * a;
  Image * b;
  Image * c;
  a = sp_image_alloc(size,size,1);
  b = sp_image_alloc(size,size,1);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = sp_image_y(a)-1; y>0;y--){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,x,y,z,sp_cinit(0,0));
	sp_image_set(b,x,y,z,sp_cinit(0,0));
      }
    }
  }  
  /*
         A      B
 
       0 1 0  0 0 0
       0 2 3  3 2 0
       0 0 0  0 1 0
   */
  /* We'll try XY reflection first */
  sp_image_set(a,1,2,0,sp_cinit(1,0));
  sp_image_set(a,1,1,0,sp_cinit(2,0));
  sp_image_set(a,2,1,0,sp_cinit(3,0));

  sp_image_set(b,1,0,0,sp_cinit(1,0));
  sp_image_set(b,1,1,0,sp_cinit(2,0));
  sp_image_set(b,0,1,0,sp_cinit(3,0));
  /* first try the out of place */
  c = sp_image_reflect(a,0,SP_AXIS_XY);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(c,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }
  /* and now inplace. Also tests reversability  */
  sp_image_reflect(b,1,SP_AXIS_XY);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }

  /* Now through the origin. Should be the same as XY for 2D images */
  sp_image_reflect(b,1,SP_ORIGO);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(c,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }
  /*
         A      B
 
       0 1 0  0 0 0
       0 2 3  0 2 3
       0 0 0  0 1 0
   */
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = sp_image_y(a)-1; y>0;y--){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(b,x,y,z,sp_cinit(0,0));
      }
    }
  }  
  sp_image_set(b,1,0,0,sp_cinit(1,0));
  sp_image_set(b,1,1,0,sp_cinit(2,0));
  sp_image_set(b,2,1,0,sp_cinit(3,0));
  sp_image_reflect(b,1,SP_AXIS_X);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }
  /*
         A      B
 
       0 1 0  0 1 0
       0 2 3  3 2 0
       0 0 0  0 0 0
   */
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = sp_image_y(a)-1; y>0;y--){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(b,x,y,z,sp_cinit(0,0));
      }
    }
  }  
  sp_image_set(b,1,2,0,sp_cinit(1,0));
  sp_image_set(b,1,1,0,sp_cinit(2,0));
  sp_image_set(b,0,1,0,sp_cinit(3,0));
  sp_image_reflect(b,1,SP_AXIS_Y);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }


}

void test_sp_image_phase_match(CuTest * tc){
  int size = 3;
  Image * a;
  Image * b;
  a = sp_image_alloc(size,size,1);
  b = sp_image_alloc(size,size,1);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = sp_image_y(a)-1; y>0;y--){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,x,y,z,sp_cinit(0,0));
	sp_image_set(b,x,y,z,sp_cinit(0,0));
      }
    }
  }  
  /*
         A      B
 
       0 1 0  0 1 0
       0 2 3  0 2 3 + pi/2 rotation
       0 0 0  0 0 0
   */
  sp_image_set(a,1,2,0,sp_cinit(1,0));
  sp_image_set(a,1,1,0,sp_cinit(2,0));
  sp_image_set(a,2,1,0,sp_cinit(3,0));

  sp_image_set(b,1,2,0,sp_cinit(0,1));
  sp_image_set(b,1,1,0,sp_cinit(0,2));
  sp_image_set(b,2,1,0,sp_cinit(0,3));
  sp_image_phase_match(a,b,1);

  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON));
      }
    }
  }

  /* Now we're gonna try with random images and a random angle perturbed with 0.1 */
  real dphi = p_drand48();
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,x,y,z,sp_cinit(p_drand48(),p_drand48()));
	sp_image_set(b,x,y,z,sp_crot(sp_image_get(a,x,y,z),dphi+p_drand48()*0.1));
	
      }
    }
  }

  sp_image_phase_match(a,b,1);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),0.1);
      }
    }
  }
  
}

void test_sp_background_adaptative_mesh(CuTest * tc){
  int size = 128;
  Complex one = sp_cinit(1,0);
  int radius = 5;
  int cx = size/2;
  int cy = size/2;
  Image * a;
  real bg_value = 10;
  a = sp_image_alloc(size,size,1);
  /* draw a circle */
  for(int x = 0;x<sp_image_x(a);x++){
    for(int y = 0;y<sp_image_y(a);y++){
      if((cx-x)*(cx-x) + (cy-y)*(cy-y) < radius*radius){
	sp_image_set(a,x,y,0,one);
      }
    }
  }
  sp_image_write(a,"background_adaptative_mesh_square.png",COLOR_GRAYSCALE);
  Image * pattern = sp_image_fft(a);
  Image * pattern_centered = sp_image_shift(pattern);
  sp_image_free(pattern);
  pattern = pattern_centered;
  pattern->scaled = 1;
  sp_image_to_intensities(pattern);
  for(int x = 0;x<sp_image_x(pattern);x++){
    for(int y = 0;y<sp_image_y(pattern);y++){
      Complex v= sp_image_get(pattern,x,y,0);
      sp_real(v) += bg_value;
      sp_image_set(pattern,x,y,0,v);
    }
  }
  sp_image_write(pattern,"background_adaptative_mesh_pattern.tif",0);
  Image * background = sp_background_adaptative_mesh(pattern,3,3,3);
  sp_image_write(background,"background_adaptative_mesh_background.tif",0);
  for(int x = 0;x<sp_image_x(pattern);x++){
    for(int y = 0;y<sp_image_y(pattern);y++){
      CuAssertDblEquals(tc,sp_cabs(sp_image_get(background,x,y,0)),bg_value,0.1);
    }
  }
}

void test_sp_image_convolute_fractional(CuTest * tc){
  int size = 10;
  Image * a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),p_drand48());
  }
  a->phased = 1;
  int precision = 1;
  Image * fcorr = sp_image_convolute_fractional(a,a,NULL,precision,1);
  Image * corr = sp_image_cross_correlate(a,a,NULL);
  for(int i = 0;i<sp_image_size(a);i++){
    CuAssertComplexEquals(tc,fcorr->image->data[i],corr->image->data[i],fabs(REAL_EPSILON*(sp_cabs(fcorr->image->data[i]))+REAL_EPSILON));  
  }
  sp_image_free(fcorr);
  Image * fconv = sp_image_convolute_fractional(a,a,NULL,precision,0);
  Image * conv = sp_image_convolute(a,a,NULL);
  for(int i = 0;i<sp_image_size(a);i++){
    CuAssertComplexEquals(tc,fconv->image->data[i],conv->image->data[i],fabs(REAL_EPSILON*(sp_cabs(fconv->image->data[i]))+REAL_EPSILON));  
  }
  sp_image_free(fconv);
  precision = 2;
  for(precision = 2;precision<10;precision++){
    fcorr = sp_image_convolute_fractional(a,a,NULL,precision,1);
    for(int z = 0;z<sp_image_z(corr);z++){
      for(int y = 0;y<sp_image_y(corr);y++){
	for(int x = 0;x<sp_image_x(corr);x++){	
	  CuAssertComplexEquals(tc,sp_image_get(fcorr,x*precision,y*precision,z*precision),sp_image_get(corr,x,y,z),1000*fabs(REAL_EPSILON*(sp_cabs(sp_image_get(fcorr,x,y,z)))+REAL_EPSILON));  
	}
      }
    }
    sp_image_free(fcorr);
    fconv = sp_image_convolute_fractional(a,a,NULL,precision,0);
    for(int z = 0;z<sp_image_z(corr);z++){
      for(int y = 0;y<sp_image_y(corr);y++){
	for(int x = 0;x<sp_image_x(corr);x++){	
	  CuAssertComplexEquals(tc,sp_image_get(fconv,x*precision,y*precision,z*precision),sp_image_get(conv,x,y,z),1000*fabs(REAL_EPSILON*(sp_cabs(sp_image_get(fconv,x,y,0)))+REAL_EPSILON));  
	}
      }
    }
    sp_image_free(fconv);
  }
  Image * b = sp_image_duplicate(a,SP_COPY_ALL);
  real t_x = 1;
  real t_y = 1;
  real t_z = 0;
  sp_image_translate(b,t_x,t_y,t_z,SP_TRANSLATE_WRAP_AROUND);
  precision = 2;
  fcorr = sp_image_convolute_fractional(a,b,NULL,precision,1);
  int x,y,z;
  long long index;
  real max = sp_image_max(fcorr,&index,&x,&y,&z);
  CuAssertTrue(tc, ((int)((real)x/precision+t_x+0.5)) % sp_image_x(a) == 0);  
  CuAssertTrue(tc, ((int)((real)y/precision+t_y+0.5)) % sp_image_y(a) == 0);  
  CuAssertTrue(tc, ((int)((real)z/precision+t_z+0.5)) % sp_image_z(a) == 0);  

  sp_image_free(b);
  b = sp_image_duplicate(a,SP_COPY_ALL);
  t_x = 0.5;
  t_y = 0;
  t_z = 0;

  sp_image_fourier_translate(b,t_x,t_y,t_z);
  precision = 2;
  fcorr = sp_image_convolute_fractional(a,b,NULL,precision,1);
  max = sp_image_max(fcorr,&index,&x,&y,&z);
  CuAssertTrue(tc, ((int)((real)x/precision+t_x+0.5)) % sp_image_x(a) == 0);  
  CuAssertTrue(tc, ((int)((real)y/precision+t_y+0.5)) % sp_image_y(a) == 0);  
  CuAssertTrue(tc, ((int)((real)z/precision+t_z+0.5)) % sp_image_z(a) == 0);  

}  


void test_sp_image_superimpose_fractional(CuTest * tc){
  /* First a simple test without enantiomorph */
  int size = 4;
  Image * a;
  Image * b;
  int precision = 1;
  a = sp_image_alloc(size,size,1);
  b = sp_image_alloc(size,size,1);
  a->phased = 1;
  b->phased = 1;
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),0);
  }
  b = sp_image_duplicate(a,SP_COPY_ALL);
  sp_image_translate(b,1,-1,0,SP_TRANSLATE_WRAP_AROUND);
  sp_image_superimpose_fractional(a,b,0,precision);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(1000*REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+1000*REAL_EPSILON));
      }
    }
  }
  /* Now test with enantiomorph*/
  sp_image_free(b);
  b = sp_image_duplicate(a,SP_COPY_ALL);
  sp_image_translate(b,1,-1,0,SP_TRANSLATE_WRAP_AROUND);
  sp_image_reflect(b,1,SP_ORIGO);
  

  sp_image_superimpose_fractional(a,b,SP_ENANTIOMORPH,precision);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON*1000));
      }
    }
  }
  /* try fractional superimpose */
  sp_image_free(b);
  b = sp_image_duplicate(a,SP_COPY_ALL);
  sp_image_fourier_translate(b,0.5,0,0);
  precision = 2;
  sp_image_superimpose_fractional(a,b,SP_ENANTIOMORPH,precision);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON*1000));
      }
    }
  }

  sp_image_free(b);
  b = sp_image_duplicate(a,SP_COPY_ALL);
  sp_image_fourier_translate(b,1.25,2.25,0);
  precision = 8;
  sp_image_superimpose_fractional(a,b,SP_ENANTIOMORPH,precision);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,x,y,z))+REAL_EPSILON*1000000));
      }
    }
  }

}

  
CuSuite* image_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();  
  SUITE_ADD_TEST(suite, test_sp_image_edge_extend);
  SUITE_ADD_TEST(suite, test_sp_bubble_sort);
  SUITE_ADD_TEST(suite, test_sp_image_median_filter);
  SUITE_ADD_TEST(suite,test_sp_image_gaussian_blur);
  SUITE_ADD_TEST(suite,test_cube_crop);
  SUITE_ADD_TEST(suite,test_sp_image_low_pass);
  SUITE_ADD_TEST(suite,test_sp_image_h5_read_write);
  SUITE_ADD_TEST(suite,test_sp_image_get_false_color);
  //  SUITE_ADD_TEST(suite,test_sp_image_noise_estimate);
  SUITE_ADD_TEST(suite,test_sp_image_dist);
  SUITE_ADD_TEST(suite,test_sp_image_max);
  SUITE_ADD_TEST(suite,test_sp_image_superimpose);
  SUITE_ADD_TEST(suite,test_sp_image_reflect);
  SUITE_ADD_TEST(suite,test_sp_image_phase_match);
  SUITE_ADD_TEST(suite,test_sp_background_adaptative_mesh);
  SUITE_ADD_TEST(suite,test_sp_image_shift);
  SUITE_ADD_TEST(suite,test_sp_image_convolute_fractional);
  SUITE_ADD_TEST(suite,test_sp_image_superimpose_fractional);
  return suite;
}

#endif
