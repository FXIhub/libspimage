#include <spimage.h>
#include <thrust/sort.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/pair.h>
#include <thrust/extrema.h>
#include <thrust/partition.h>
#include <thrust/count.h>
#include <cstdlib>

static real bezier_map_interpolation(sp_smap * map, real x);
static void support_from_absolute_threshold_cuda(SpPhaser * ph, cufftComplex * blur, real abs_threshold);

static __global__ void CUDA_support_from_threshold(cufftComplex * a,float abs_threshold,int * pixel_flags, int size){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    if(a[i].x > abs_threshold){
      pixel_flags[i] |= SpPixelInsideSupport;
    }else{
      pixel_flags[i] &= ~SpPixelInsideSupport;
    }
  }
}


static __global__ void CUDA_dephase(cufftComplex * a, int size){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
#ifdef _STRICT_IEEE_754
    a[i].x = __fsqrt_rn(__fadd_rn(__fmul_rn(a[i].x,a[i].x), __fmul_rn(a[i].y,a[i].y)));
#else
    a[i].x = sqrt(a[i].x*a[i].x+a[i].y*a[i].y);
#endif
    a[i].y = 0;
  }
}

static __global__ void CUDA_complex_to_float(cufftComplex * in,float * out, int size){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    out[i] = in[i].x;
  }
}

struct maxCufftComplex{   
  __host__ __device__ cufftComplex operator()(const cufftComplex lhs, const cufftComplex rhs) { 
    if(lhs.x*lhs.x + lhs.y*lhs.y < rhs.x*rhs.x + rhs.y*rhs.y){
      return rhs;
    }
    return lhs;
  }
};

struct descend_cmpCufftComplex{   
   __device__ bool operator()(const cufftComplex lhs, const cufftComplex rhs) { 
#ifdef _STRICT_IEEE_754
    if(__fadd_rn(__fmul_rn(lhs.x,lhs.x), __fmul_rn(lhs.y,lhs.y))
       > __fadd_rn(__fmul_rn(rhs.x,rhs.x), __fmul_rn(rhs.y,rhs.y))){
      return true;
    }
#else
    if(lhs.x*lhs.x + lhs.y*lhs.y > rhs.x*rhs.x + rhs.y*rhs.y){
      return true;
    }
#endif
    return false;
  }
};

struct descend_real_cmpCufftComplex{   
   __device__ bool operator()(const cufftComplex lhs, const cufftComplex rhs) const { 
#ifdef _STRICT_IEEE_754
    if(__fadd_rn(__fmul_rn(lhs.x,lhs.x), __fmul_rn(lhs.y,lhs.y))
       > __fadd_rn(__fmul_rn(rhs.x,rhs.x), __fmul_rn(rhs.y,rhs.y))){
      return true;
    }
#else
    if(lhs.x > rhs.x){
      return true;
    }
#endif
    return false;
  }
};


struct is_larger_than
{
  float t;
  is_larger_than(float _t)
  :t(_t)
  {
  }
   __device__  bool operator()(const float x) const
  {
    return (x > t);
  }
};

struct is_larger_than_Complex
{
  float t;
  is_larger_than_Complex(cufftComplex _t)
  :t(_t.x)
  {
  }
   __device__  bool operator()(const cufftComplex x) const
  {
    return (x.x > t);
  }
};
/*
static __global__ void CUDA_simple_histogram(float * a, int size,float * bin_ceil, int nbins, int * output){
  __shared__ int data[512];
  // initialize shared storage
  data[threadIdx.x] = 0;
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<size){
    for(int j = 0;j<nbins;j++){
      if(bin_ceil[j] > a[i]){
	data[j]++;
      }
    }
  }
  __syncthreads();
  if(threadIdx.x < nbins){
    output[threadIdx.x] += data[threadIdx.x];
  }
}

static float sp_get_highest_n_cuda3(float * a, int n, int size){
  thrust::device_ptr<float> begin =  thrust::device_pointer_cast(a);
  thrust::device_ptr<float> end =  thrust::device_pointer_cast(a+size);
  thrust::device_ptr<float> bottom = begin;
  thrust::device_ptr<float> top = end;
  thrust::pair<thrust::device_ptr<float>,thrust::device_ptr<float> > min_max = thrust::minmax_element(bottom,top);  
  float min = *(min_max.first);
  float max = *(min_max.second);  
  int bins = 256;
  float * bin_ceil = (float *)malloc(sizeof(float)*bins);
  for(int i = 0;i<bins;i++){
    bin_ceil[i] = min+(max-min)*(i+1.0f)/bins;
  }
  float * d_bin_ceil;
  cutilSafeCall(cudaMalloc((void **)&d_bin_ceil,sizeof(float)*bins));
  cutilSafeCall(cudaMemcpy(d_bin_ceil,bin_ceil,sizeof(float)*bins,cudaMemcpyHostToDevice));
  int * output = (int *)malloc(sizeof(int)*bins);
  int * d_output;
  cutilSafeCall(cudaMalloc((void **)&d_output,sizeof(int)*bins));
  int threads_per_block = 512;
  int number_of_blocks = (size+threads_per_block-1)/threads_per_block;  
    CUDA_simple_histogram<<<number_of_blocks, threads_per_block>>>(a,size,d_bin_ceil,bins,d_output);
  cutilSafeCall(cudaMemcpy(output,d_output,sizeof(int)*bins,cudaMemcpyDeviceToHost));
}

static float sp_get_highest_n_cuda(cufftComplex * a, int n, int size){
  thrust::device_ptr<cufftComplex> begin =  thrust::device_pointer_cast(a);
  thrust::device_ptr<cufftComplex> end =  thrust::device_pointer_cast((cufftComplex *)(a+size));
  thrust::device_ptr<cufftComplex> bottom = begin;
  thrust::device_ptr<cufftComplex> top = end;
  thrust::pair<thrust::device_ptr<cufftComplex>,thrust::device_ptr<cufftComplex> > min_max = thrust::minmax_element(bottom,top,descend_real_cmpCufftComplex());  
  thrust::device_ptr<cufftComplex> min_p = thrust::min_element(begin,end,descend_real_cmpCufftComplex());  
  cufftComplex min = *(min_max.first);
  cufftComplex max = *(min_max.second);  
  cufftComplex guess = {0,0};
  do{
    guess.x= (min.x+max.x)/2.0f;
    int my_n = thrust::count_if(begin,end,is_larger_than_Complex(guess));
    if(my_n < n){
      if(max.x == guess.x){
	return guess.x;
      }
      max = guess;
    }else if(my_n > n){
      if(min.x == guess.x){
	return guess.x;
      }
      min = guess;
    }else{
      return guess.x;
    }      
  }while(1);  
}


static float sp_get_highest_n_cuda4(float * a, int n, int size){
  thrust::device_ptr<float> begin =  thrust::device_pointer_cast(a);
  thrust::device_ptr<float> end =  thrust::device_pointer_cast(a+size);
  //  thrust::device_ptr<float> bottom = begin;
  //  thrust::device_ptr<float> top = end;
  thrust::pair<thrust::device_ptr<float>,thrust::device_ptr<float> > min_max = thrust::minmax_element(begin,end);  
  float min = *(min_max.first);
  float max = *(min_max.second);  
  int my_n = size/2;
  int bottom = 0;
  int top = size;
  do{
    float guess = (min+max)/2.0f;
    my_n = thrust::count_if(begin,end,is_larger_than(guess));
    if(my_n < n){
      if(max == guess){
	return guess;
      }
      max = guess;
    }else if(my_n > n){
      if(min == guess){
	return guess;
      }
      min = guess;
    }else{
      return guess;
    }      
  }while(1);  
}


static float sp_get_highest_n_cuda2(float * a, int n, int size){
  thrust::device_ptr<float> begin =  thrust::device_pointer_cast(a);
  thrust::device_ptr<float> end =  thrust::device_pointer_cast(a+size);
  thrust::device_ptr<float> bottom = begin;
  thrust::device_ptr<float> top = end;
  while(top-bottom > 1024*16){
    thrust::pair<thrust::device_ptr<float>,thrust::device_ptr<float> > min_max = thrust::minmax_element(bottom,top);  
    float min = *(min_max.first);
    float max = *(min_max.second);
    float guess = (min+max)/2;
    thrust::device_ptr<float> middle = thrust::partition(bottom,top,is_larger_than(guess));
    int first =  middle-begin;
    if(first < n){
      bottom = middle;
    }else if(first > n){
      top = middle;
    }else{
      return *middle;
    }      
  }  
  thrust::sort(bottom,top);
  return begin[n];
}
*/

static float sp_image_max_cuda(cufftComplex * a, int size){
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(a);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(a+size));
  cufftComplex max = {0,0};
  max = thrust::reduce(beginc,endc,max,maxCufftComplex());
  float ret = max.x*max.x+max.y*max.y;
  return sqrt(ret);
}
/*
static int descend_complex_compare(const void * pa,const void * pb){
  Complex a,b;
  a = *((Complex *)pa);
  b = *((Complex *)pb);
  if(sp_cabs(a) < sp_cabs(b)){
    return 1;
  }else if(sp_cabs(a) == sp_cabs(b)){
    return 0;
  }else{
    return -1;
  }
}
*/
/*
static void sp_image_sort_and_test_cuda(cufftComplex * a, int size){
  Image * test_golden = sp_image_alloc(size,1,1);
  cutilSafeCall(cudaMemcpy(test_golden->image->data,a,sizeof(cufftComplex)*size,cudaMemcpyDeviceToHost));
  qsort(test_golden->image->data,size,sizeof(Complex),descend_complex_compare);
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(a);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(a+size));
  thrust::sort(beginc,endc,descend_cmpCufftComplex());
  Image * test_gpu = sp_image_alloc(size,1,1);
  cutilSafeCall(cudaMemcpy(test_gpu->image->data,a,sizeof(cufftComplex)*size,cudaMemcpyDeviceToHost));
  for(int i =0 ;i<size;i++){
    if(sp_cabs(test_golden->image->data[i]) != sp_cabs(test_gpu->image->data[i])){
      fprintf(stderr,"error!\n");
    }
  }
}

static void sp_image_sort_cuda(cufftComplex * a, int size){
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(a);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(a+size));
  thrust::sort(beginc,endc,descend_real_cmpCufftComplex());
}
*/

int sp_support_area_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph){
  SpSupportAreaParameters * params = (SpSupportAreaParameters *)alg->params;
  real radius =  bezier_map_interpolation(params->blur_radius_map,ph->iteration);
  cufftComplex * blur;
  cutilSafeCall(cudaMalloc((void**)&blur, sizeof(cufftComplex)*ph->image_size));
  float * sort;
    cutilSafeCall(cudaMalloc((void**)&sort, sizeof(float)*ph->image_size));
  cutilSafeCall(cudaMemcpy(blur,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToDevice));
  /* we need to dephase the image first */
  CUDA_dephase<<<ph->number_of_blocks, ph->threads_per_block>>>(blur,ph->image_size);
  sp_cuda_check_errors();
  sp_gaussian_blur_cuda(blur,blur,sp_3matrix_x(ph->amplitudes),sp_3matrix_y(ph->amplitudes),sp_3matrix_z(ph->amplitudes),radius,ph->cufft_plan);

  CUDA_complex_to_float<<<ph->number_of_blocks, ph->threads_per_block>>>(blur,sort,ph->image_size);
  //  cutilSafeCall(cudaMemcpy(sort,blur,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToDevice));
  real area = bezier_map_interpolation(params->area,ph->iteration);

  thrust::device_ptr<float> begin =  thrust::device_pointer_cast(sort);
  thrust::device_ptr<float> end =  thrust::device_pointer_cast(sort+ph->image_size);
  
  thrust::sort(begin,end);
  /* sort in ascending order, that's the reason of 1.0-area */
  float v = 0;
  cutilSafeCall(cudaMemcpy(&v,&sort[(int)(ph->image_size*(1.0-area))],sizeof(float),cudaMemcpyDeviceToHost));
  
  real abs_threshold = v;
  support_from_absolute_threshold_cuda(ph,blur,abs_threshold);
  cutilSafeCall(cudaFree(blur));
  cutilSafeCall(cudaFree(sort));
  return 0;
}

int sp_support_template_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser *ph){
  SpSupportTemplateParameters *params = (SpSupportTemplateParameters *)alg->params;
  cufftComplex *blured;
  cutilSafeCall(cudaMalloc((void**)&blured, sizeof(cufftComplex)*ph->image_size));
  cutilSafeCall(cudaMemcpy(blured,params->blured->image->data,sizeof(cufftComplex)*ph->image_size,cudaMemcpyHostToDevice)); //should this be device to host?
  real area = bezier_map_interpolation(params->area,ph->iteration);
  real abs_threshold = sp_cabs(params->sorted->image->data[(int)(sp_image_size(params->sorted)*area*params->original_area)]);
  support_from_absolute_threshold_cuda(ph,blured,abs_threshold);
  cutilSafeCall(cudaFree(blured));
  return 0;
}

int sp_support_threshold_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph){
  SpSupportThresholdParameters * params = (SpSupportThresholdParameters *)alg->params;
  real radius =  bezier_map_interpolation(params->blur_radius_map,ph->iteration);
  cufftComplex * blur;
  cutilSafeCall(cudaMalloc((void**)&blur, sizeof(cufftComplex)*ph->image_size));
  cutilSafeCall(cudaMemcpy(blur,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToDevice));
  /* we need to dephase the image first */
  CUDA_dephase<<<ph->number_of_blocks, ph->threads_per_block>>>(blur,ph->image_size);
  sp_cuda_check_errors();
  sp_gaussian_blur_cuda(blur,blur,sp_3matrix_x(ph->amplitudes),sp_3matrix_y(ph->amplitudes),sp_3matrix_z(ph->amplitudes),radius,ph->cufft_plan);
  real rel_threshold = bezier_map_interpolation(params->threshold,ph->iteration);  
  real abs_threshold = sp_image_max_cuda(blur,ph->image_size)*rel_threshold;  
  support_from_absolute_threshold_cuda(ph,blur,abs_threshold);
  cutilSafeCall(cudaFree(blur));
  return 0;
}

int sp_support_static_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph){
  return 0;
}

int sp_support_close_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph){
  SpSupportCloseParameters *params = alg->params;
  cutilSafeCall(cudaMemcpy(pixel_flags->data,d_pixel_flags,sizeof(int)*ph->image_size,cudaMemcpyDeviceToHost));
  sp_i3matrix *in = ph->pixel_flags;
  sp_i3matrix *tmp1 = sp_i3matrix_duplicate(ph->pixel_flags);
  sp_i3matrix *tmp2 = sp_i3matrix_duplicate(ph->pixel_flags);
  sp_i3matrix *foo;

  for (int i = 0; i < params->size; i++) {
    for (int x = 0; x < sp_i3matrix_x(in); x++) {
      for (int y = 0; y < sp_i3matrix_y(in); y++) {
	for (int z = 0; z < sp_i3matrix_z(in); z++) {
	  if ((sp_i3matrix_get(tmp1,x,y,z) == SpPixelInsideSupport) ||
	      ((x != 0 && sp_i3matrix_get(tmp1,x-1,y,z) ==
		SpPixelInsideSupport) ||
	       (x != sp_image_x(tmp1)-1 && sp_i3matrix_get(tmp1,x+1,y,z) ==
		SpPixelInsideSupport) ||
	       (y != 0 && sp_i3matrix_get(tmp1,x,y-1,z) == 
		SpPixelInsideSupport) ||
	       (y != sp_image_y(tmp1)-1 && sp_i3matrix_get(tmp1,x,y+1,z) ==
		SpPixelInsideSupport) ||
	       (z != 0 && sp_i3matrix_get(tmp1,x,y,z-1) ==
		SpPixelInsideSupport) ||
	       (z != sp_image_z(tmp1)-1 && sp_i3matrix_get(tmp1,x,y,z+1) ==
		SpPixelInsideSupport))) {
	    tmp2->data[z*sp_i3matrix_y(in)*sp_i3matrix_x(in)+
		       y*sp_i3matrix_x(in)+x] |= SpPixelInsideSupport;
	  } else {
	    tmp2->data[z*sp_i3matrix_y(in)*sp_i3matrix_x(in)+
		       y*sp_i3matrix_x(in)+x] &= ~SpPixelInsideSupport;
	  }
	}
      }
    }
    foo = tmp1;
    tmp1 = tmp2;
    tmp2 = foo;
  }

  for (int i = 0; i < params->size; i++) {
    for (int x = 0; x < sp_i3matrix_x(in); x++) {
      for (int y = 0; y < sp_i3matrix_y(in); y++) {
	for (int z = 0; z < sp_i3matrix_z(in); z++) {
	  if ((sp_i3matrix_get(tmp1,x,y,z) == ~SpPixelInsideSupport) ||
	      ((x != 0 && sp_i3matrix_get(tmp1,x-1,y,z) ==
		~SpPixelInsideSupport) ||
	       (x != sp_image_x(tmp1)-1 && sp_i3matrix_get(tmp1,x+1,y,z) ==
		~SpPixelInsideSupport) ||
	       (y != 0 && sp_i3matrix_get(tmp1,x,y-1,z) == 
		~SpPixelInsideSupport) ||
	       (y != sp_image_y(tmp1)-1 && sp_i3matrix_get(tmp1,x,y+1,z) ==
		~SpPixelInsideSupport) ||
	       (z != 0 && sp_i3matrix_get(tmp1,x,y,z-1) ==
		~SpPixelInsideSupport) ||
	       (z != sp_image_z(tmp1)-1 && sp_i3matrix_get(tmp1,x,y,z+1) ==
		~SpPixelInsideSupport))) {
	    tmp2->data[z*sp_i3matrix_y(in)*sp_i3matrix_x(in)+
		       y*sp_i3matrix_x(in)+x] &= ~SpPixelInsideSupport;
	  } else {
	    tmp2->data[z*sp_i3matrix_y(in)*sp_i3matrix_x(in)+
		       y*sp_i3matrix_x(in)+x] |= SpPixelInsideSupport;
	  }
	}
      }
    }
    foo = tmp1;
    tmp1 = tmp2;
    tmp2 = foo;
  }
  sp_i3matrix_free(tmp1);
  sp_i3matrix_free(tmp2);

  return 0;
  cutilSafeCall(cudaMemcpy(d_pixel_flags,d_pixel_flags->data,sizeof(int)*ph->image_size,cudaMemcpyHostToDevice));
}

static void support_from_absolute_threshold_cuda(SpPhaser * ph, cufftComplex * blur, real abs_threshold){
  CUDA_support_from_threshold<<<ph->number_of_blocks, ph->threads_per_block>>>(blur,abs_threshold,ph->d_pixel_flags,ph->image_size);
  sp_cuda_check_errors();
}

static real bezier_map_interpolation(sp_smap * map, real x){
  sp_list * keys = sp_smap_get_keys(map);
  sp_list * values = sp_smap_get_values(map);
  unsigned int idx = 0;
  for(idx = 0;idx<sp_list_size(keys);idx++){
    if(x < sp_list_get(keys,idx)){
      break;
    }
  }
  if(idx == 0){
    return sp_list_get(values,0);
  }
  if(idx == sp_list_size(keys)){
    return sp_list_get(values,sp_list_size(keys)-1);
  }
  /* Cubic Bezier curve taken from http://en.wikipedia.org/wiki/Bézier_curve */
  real p0y = sp_list_get(values,idx-1);
  real p0x = sp_list_get(keys,idx-1);
  real p3y = sp_list_get(values,idx);
  real p3x = sp_list_get(keys,idx);
  real t = (2 - 2*(p0x - p3x)*
	    pow(8*p0x*p3x*x - 2*p3x*pow(p0x,2) - 4*x*pow(p0x,2) + 
		2*pow(p0x,3) - 2*p0x*pow(p3x,2) - 4*x*pow(p3x,2) + 
		2*pow(p3x,3) + pow(pow(p0x - p3x,4)*
				   (6*p0x*p3x - 16*p0x*x - 16*p3x*x + 5*pow(p0x,2) +
				    5*pow(p3x,2) + 16*pow(x,2)),0.5f),
		-0.3333333333333333f) +
	    2*pow(p0x - p3x,-1)*
	    pow(8*p0x*p3x*x - 2*p3x*pow(p0x,2) - 4*x*pow(p0x,2) +
		2*pow(p0x,3) - 2*p0x*pow(p3x,2) - 4*x*pow(p3x,2) +
		2*pow(p3x,3) + pow(pow(p0x - p3x,4)*
				   (6*p0x*p3x - 16*p0x*x -
				    16*p3x*x + 5*pow(p0x,2) + 
				     5*pow(p3x,2) + 16*pow(x,2)),0.5f)
		,0.3333333333333333f)
	    )/4.;  
  return 3*p0y*t*pow(1 - t,2) + p0y*pow(1 - t,3) +
    3*p3y*(1 - t)*pow(t,2) + p3y*pow(t,3);
}
