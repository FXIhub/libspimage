#ifndef PNG_DEBUG
#  define PNG_DEBUG 3
#endif

#define _XOPEN_SOURCE 500

#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <tiffio.h>
#include <png.h>
#include <float.h>
#include <ctype.h>
#include <strings.h>
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif

#include "spimage.h"

static void write_h5_img(const Image * img,const char * filename, int output_precision);
static Image * _read_imagefile(const char * filename,const char * file, int line);
static Image * read_tiff(const char * filename);
static  void write_tiff(const Image * img,const char * filename);
static  void write_csv(const Image * img,const char * filename);
static Image * read_png(const char * filename);
static int write_png(const Image * img,const char * filename, int color);
static int write_vtk(const Image * img,const char * filename);
static int write_xplor(const Image * img,const char * filename);
static Image * read_smv(const char * filename);
static Image * read_anton_datafile(hid_t file_id,hid_t dataset_id, const char * filename);
static void write_cxdi(const Image * img,const char * filename);
static Image * read_cxdi(const char * filename);

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


void sp_image_write(const Image * img, const const char * filename, int flags){
  char buffer[1024];
  strcpy(buffer,filename);
  for(int i = 0;i<strlen(buffer);i++){
    buffer[i] = tolower(buffer[i]);
  }
 /* select the correct function depending on the buffer extension */
  if(strrchr(buffer,'.') && strcmp(strrchr(buffer,'.'),".h5") == 0){
    /* we have an h5 file */
    write_h5_img(img,filename,sizeof(real));
  }else if(strrchr(buffer,'.') && strcmp(strrchr(buffer,'.'),".png") == 0){
    write_png(img,filename,flags);
  }else if(strrchr(buffer,'.') && strcmp(strrchr(buffer,'.'),".vtk") == 0){
    write_vtk(img,filename);
  }else if(strrchr(buffer,'.') && (strcmp(strrchr(buffer,'.'),".tif") == 0 ||strcmp(strrchr(buffer,'.'),".tiff") == 0 )){
    write_tiff(img,filename);
  }else if(strrchr(buffer,'.') && (strcmp(strrchr(buffer,'.'),".csv") == 0)){
    if(img->num_dimensions == SP_3D){
      fprintf(stderr,"Cannot export 3D file to csv");
    }
    write_csv(img,filename);
  }else if(strrchr(buffer,'.') && (strcmp(strrchr(buffer,'.'),".xplor") == 0)){
    if(img->num_dimensions != SP_3D){
      fprintf(stderr,"Can only export 3D files to xplor");
    }
    write_xplor(img,filename);
  }else if(strrchr(buffer,'.') && (strcmp(strrchr(buffer,'.'),".cxdi") == 0)){
    write_cxdi(img,filename);
  }else{
    fprintf(stderr,"Unsupported file type: %s\n",filename);
    abort();
  }
}

Image * _sp_image_read(const char * filename, int flags, const char * file, int line){
  char buffer[1024];
  strcpy(buffer,filename);
  for(int i = 0;i<strlen(buffer);i++){
    buffer[i] = tolower(buffer[i]);
  }
  /* select the correct function depending on the filename extension */
  if(strrchr(buffer,'.') && strcmp(strrchr(buffer,'.'),".h5") == 0){
    /* we have an h5 file */
    return _read_imagefile(filename,file,line);
  }else if(strrchr(buffer,'.') && strcmp(strrchr(buffer,'.'),".png") == 0){
    /* we  have a png file */
    return read_png(filename);
  }else if(strrchr(buffer,'.') && strcmp(strrchr(buffer,'.'),".vtk") == 0){
    /* we have a vtk file */
    fprintf(stderr,"Cannot read VTK files!\n");
    return NULL;
  }else if(strrchr(buffer,'.') && (strcmp(strrchr(buffer,'.'),".tif") == 0 ||strcmp(strrchr(buffer,'.'),".tiff") == 0 )){
    /* we have a tiff file */
    return read_tiff(filename);
  }else if(strrchr(buffer,'.') && (strcmp(strrchr(buffer,'.'),".smv") == 0 ||strcmp(strrchr(buffer,'.'),".SMV") == 0 )){
    /* we have an smv file */
    return read_smv(filename);
  }else if(strrchr(buffer,'.') && (strcmp(strrchr(buffer,'.'),".cxdi") == 0 ||strcmp(strrchr(buffer,'.'),".CXDI") == 0 )){
    /* we have an hdf5 simple data file file */
    return read_cxdi(filename);
  }else{
    fprintf(stderr,"Unsupported file type: %s\n",filename);
    abort();
  }
  return NULL;
}


/* superseeded by the new one. Just leave it here in case strange bugs show up in the new one */
/*
static void write_h5_img(Image * img,const char * filename, int output_precision){
  hid_t dataspace_id;
  hid_t dataset_id;
  hid_t file_id;
  int status;
  int version;
  hsize_t  dims[2];
  real values[2];
  sp_3matrix * tmp;
  int i;
  hid_t out_type_id = 0;
  hid_t mem_type_id = 0;
  hid_t plist;
  hsize_t chunk_size[2] = {sp_c3matrix_x(img->image),sp_c3matrix_y(img->image)};
  if(output_precision == sizeof(double)){
    out_type_id = H5T_NATIVE_DOUBLE;
  }else if(output_precision == sizeof(float)){
    out_type_id = H5T_NATIVE_FLOAT;
  }else{
    abort();
  }
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }

  dims[0] = sp_c3matrix_x(img->image);
  dims[1] = sp_c3matrix_y(img->image);
  file_id = H5Fcreate(filename,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple( 2, dims, NULL );

  plist = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk(plist,2,chunk_size);
  H5Pset_deflate(plist,6);

  dataset_id = H5Dcreate(file_id, "/mask", H5T_NATIVE_INT,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id,H5T_NATIVE_INT , H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, img->mask->data);
  status = H5Dclose(dataset_id);

  tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),1);
  for(i = 0;i<sp_image_size(img);i++){
    tmp->data[i] = sp_real(img->image->data[i]);
  }

  dataset_id = H5Dcreate(file_id, "/real", out_type_id,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
  status = H5Dclose(dataset_id);
  sp_3matrix_free(tmp);

  if(img->phased){
    tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),1);
    for(i = 0;i<sp_image_size(img);i++){
      tmp->data[i] = sp_imag(img->image->data[i]);
    }

    dataset_id = H5Dcreate(file_id, "/imag",out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
    status = H5Dclose(dataset_id);
    sp_3matrix_free(tmp);

  }
  dims[0] = 2;
  dataspace_id = H5Screate_simple( 1, dims, NULL );
  dataset_id = H5Dcreate(file_id, "/image_center",out_type_id ,
			 dataspace_id, H5P_DEFAULT);
  values[0] = img->detector->image_center[0];
  values[1] = img->detector->image_center[1];
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

  dims[0] = 1;
  dataspace_id = H5Screate_simple( 1, dims, NULL );
  values[0] = img->phased;
  dataset_id = H5Dcreate(file_id, "/phased", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->shifted;
  dataset_id = H5Dcreate(file_id, "/shifted", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->detector->wavelength;
  dataset_id = H5Dcreate(file_id, "/wavelength", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->detector->pixel_size[0];
  values[1] = img->detector->pixel_size[1];
  values[2] = img->detector->pixel_size[2];
  dataset_id = H5Dcreate(file_id, "/pixel_size", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->detector->detector_distance;
  dataset_id = H5Dcreate(file_id, "/detector_distance", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->scaled;
  dataset_id = H5Dcreate(file_id, "/scaled", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);


  version = 2;
  dataset_id = H5Dcreate(file_id, "/version", H5T_NATIVE_INT,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, &version);
  status = H5Dclose(dataset_id);


  status = H5Sclose(dataspace_id);


  status = H5Fclose(file_id);
}
*/

static void write_h5_img(const Image * img,const char * filename, int output_precision){
  hid_t dataspace_id;
  hid_t dataset_id;
  hid_t file_id;
  int status;
  int version;
  hsize_t  dims[3];
  real values[3];
  sp_3matrix * tmp;
  int i;
  hid_t out_type_id = 0;
  hid_t mem_type_id = 0;
  hid_t plist;
  hsize_t chunk_size[3] = {sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),sp_c3matrix_z(img->image)};
  char tmpfile[1024];
  //  strcpy(tmpfile,filename);
  //  strcat(tmpfile,"XXXXXX");
  sprintf(tmpfile,"%s-%d",filename,rand());
  /*  int fd = mkstemp(tmpfile);
  if(fd == -1){
    sp_error_warning("Unable create temporary filename");
    return;
  }
  close(fd);*/
  if(output_precision == sizeof(double)){
    out_type_id = H5T_NATIVE_DOUBLE;
  }else if(output_precision == sizeof(float)){
    out_type_id = H5T_NATIVE_FLOAT;
  }else{
    abort();
  }
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }

  dims[0] = sp_c3matrix_x(img->image);
  dims[1] = sp_c3matrix_y(img->image);
  dims[2] = sp_c3matrix_z(img->image);
  //  file_id = H5Fcreate(filename,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  file_id = H5Fcreate(tmpfile,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple( 3, dims, NULL );

  plist = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk(plist,3,chunk_size);
  H5Pset_deflate(plist,6);

  dataset_id = H5Dcreate(file_id, "/mask", H5T_NATIVE_INT,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id,H5T_NATIVE_INT , H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, img->mask->data);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }

  tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),
			 sp_c3matrix_z(img->image));
  for(i = 0;i<sp_image_size(img);i++){
    tmp->data[i] = sp_real(img->image->data[i]);
  }

  dataset_id = H5Dcreate(file_id, "/real", out_type_id,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }
  sp_3matrix_free(tmp);

  if(img->phased){
    tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),
			   sp_c3matrix_z(img->image));
    for(i = 0;i<sp_image_size(img);i++){
      tmp->data[i] = sp_imag(img->image->data[i]);
    }

    dataset_id = H5Dcreate(file_id, "/imag",out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
    if(status < 0){
      goto error;
    }
    status = H5Dclose(dataset_id);
    if(status < 0){
      goto error;
    }
    sp_3matrix_free(tmp);

  }
  dims[0] = 3;
  dataspace_id = H5Screate_simple( 1, dims, NULL );
  if(dataspace_id < 0){
    goto error;
  }
  dataset_id = H5Dcreate(file_id, "/image_center",out_type_id ,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  values[0] = img->detector->image_center[0];
  values[1] = img->detector->image_center[1];
  values[2] = img->detector->image_center[2];
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }
  status = H5Sclose(dataspace_id);
  if(status < 0){
    goto error;
  }

  dims[0] = 1;
  dataspace_id = H5Screate_simple( 1, dims, NULL );
  if(dataspace_id < 0){
    goto error;
  }

  values[0] = img->phased;
  dataset_id = H5Dcreate(file_id, "/phased", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }

  values[0] = img->shifted;
  dataset_id = H5Dcreate(file_id, "/shifted", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }

  values[0] = img->detector->wavelength;
  dataset_id = H5Dcreate(file_id, "/lambda", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }

  values[0] = img->detector->pixel_size[0];
  values[1] = img->detector->pixel_size[1];
  values[2] = img->detector->pixel_size[2];
  dataset_id = H5Dcreate(file_id, "/pixel_size", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }

  values[0] = img->num_dimensions;
  dataset_id = H5Dcreate(file_id, "/num_dimensions", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }

  values[0] = img->detector->detector_distance;
  dataset_id = H5Dcreate(file_id, "/detector_distance", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }

  values[0] = img->scaled;
  dataset_id = H5Dcreate(file_id, "/scaled", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }



  version = 2;
  dataset_id = H5Dcreate(file_id, "/version", H5T_NATIVE_INT,
			 dataspace_id, H5P_DEFAULT);
  if(dataset_id < 0){
    goto error;
  }
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, &version);
  if(status < 0){
    goto error;
  }
  status = H5Dclose(dataset_id);
  if(status < 0){
    goto error;
  }


  status = H5Sclose(dataspace_id);
  if(status < 0){
    goto error;
  }


  status = H5Fclose(file_id);
  if(status < 0){
    goto error;
  }
  if(rename(tmpfile,filename)){
    sp_error_warning("Unable to rename %s to %s",tmpfile,filename);
  }
  return;
  
 error:
  sp_error_warning("Error while writing HDF5 file at %s:%d\n",__FILE__,__LINE__);
  return;    
}


Image * _read_imagefile(const char * filename,const char * file, int line){
  Image * res = sp_malloc(sizeof(Image));
  hid_t file_id,dataset_id,space;
  int status,i;
  int version;
  hsize_t dims[3];
  hid_t mem_type_id = 0;
  H5E_auto_t func;
  void * client_data;
  real values[3] = {0,0,0};
  sp_3matrix * tmp;
  int flag_num_dimensions = 0;
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }
  
  
  
  res->detector = sp_malloc(sizeof(Detector));
  
  
  H5Eget_auto(&func,&client_data);
  /* turn off warning to check file and version because they might not exist */
  H5Eset_auto(NULL,NULL);  

  file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0){
    sp_error_warning("Unable to open %s",filename);
    return NULL;
  }
  

  dataset_id = H5Dopen(file_id, "/data/data");
  if(dataset_id >= 0){
    /* we have Anton's simple data format */
    H5Eset_auto(func,client_data);
    return read_anton_datafile(file_id,dataset_id,filename);
  }

  dataset_id = H5Dopen(file_id, "/version");
  /* File includes version information */
  if(dataset_id>=0){
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, &version);
    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }
    status = H5Dclose(dataset_id);
    if(version == 2){
      dataset_id = H5Dopen(file_id, "/mask");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }
      space = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(space,dims,NULL);
      if(H5Sget_simple_extent_ndims(space) == 3){
	res->image = _sp_c3matrix_alloc(dims[0],dims[1],dims[2],file,line);
	res->mask = _sp_i3matrix_alloc(dims[0],dims[1],dims[2],file,line);
      }else{
	res->image = _sp_c3matrix_alloc(dims[0],dims[1],1,file,line);
	res->mask = _sp_i3matrix_alloc(dims[0],dims[1],1,file,line);
      }
      
      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, res->mask->data);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }


      status = H5Dclose(dataset_id);
      
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      dataset_id = H5Dopen(file_id, "/image_center");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }
      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }
      status = H5Dclose(dataset_id);
      res->detector->image_center[0] = values[0];
      res->detector->image_center[1] = values[1];
      if(values[2]){
	res->detector->image_center[2] = values[2];
      }else{
	res->detector->image_center[2] = 0;
      }
      
      dataset_id = H5Dopen(file_id, "/phased");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      res->phased = values[0];
      
      dataset_id = H5Dopen(file_id, "/shifted");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);

      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      res->shifted = values[0];
      
      dataset_id = H5Dopen(file_id, "/scaled");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      res->scaled = values[0];
      
      dataset_id = H5Dopen(file_id, "/detector_distance");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id,  mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      res->detector->detector_distance = values[0];
      
      dataset_id = H5Dopen(file_id, "/lambda");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      res->detector->wavelength = values[0];
      
      dataset_id = H5Dopen(file_id, "/pixel_size");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      res->detector->pixel_size[0] = values[0];
      res->detector->pixel_size[1] = values[1];
      res->detector->pixel_size[2] = values[2];
      

      H5Eget_auto(&func,&client_data);
      /* turn off warning to check num_dimensions because it might not exist */
      H5Eset_auto(NULL,NULL);
      dataset_id = H5Dopen(file_id, "/num_dimensions");
      H5Eset_auto(func,client_data);
      if(dataset_id>=0){
	flag_num_dimensions = 1;
	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, values);
	if(status < 0){
	  sp_error_warning("Unable to read dataset from file %s",filename);
	  return NULL;
	}

	status = H5Dclose(dataset_id);
	res->num_dimensions = values[0];
      }else{
	/* we'll try to guess the dimensions */
	res->num_dimensions = 0;
	if(sp_image_x(res) > 1){
	  res->num_dimensions++;
	}
	if(sp_image_y(res) > 1){
	  res->num_dimensions++;
	}
	if(sp_image_z(res) > 1){
	  res->num_dimensions++;
	}	  
      }
      
      if(res->phased){
	tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),
			       sp_i3matrix_y(res->mask),
			       sp_i3matrix_z(res->mask),file,line);
	dataset_id = H5Dopen(file_id, "/imag");
	if(dataset_id < 0){
	  sp_error_warning("Unable to open dataset in file %s",filename);
	  return NULL;
	}

	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
	if(status < 0){
	  sp_error_warning("Unable to read dataset from file %s",filename);
	  return NULL;
	}

	status = H5Dclose(dataset_id);
	if(status < 0){
	  sp_error_warning("Unable to close dataset from file %s",filename);
	  return NULL;
	}

	for(i = 0;i<sp_3matrix_size(tmp);i++){
	  sp_imag(res->image->data[i]) = tmp->data[i];
	  sp_real(res->image->data[i]) = 0;
	}
	sp_3matrix_free(tmp);
      }
      
      tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			     sp_i3matrix_z(res->mask),file,line);
      if(!tmp){
	sp_error_warning("Unable to allocate matrix");
	return NULL;
      }

      dataset_id = H5Dopen(file_id, "/real");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      for(i = 0;i<sp_3matrix_size(tmp);i++){
	sp_real(res->image->data[i]) += tmp->data[i];
      }
      sp_3matrix_free(tmp);
      
      status = H5Fclose(file_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }
    }
  }else{
    /* File does *NOT* includes version information */
    dataset_id = H5Dopen(file_id, "/mask");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }

    space = H5Dget_space(dataset_id);
    if(space < 0){
      sp_error_warning("Unable to get space in file %s",filename);
      return NULL;
    }
    res->num_dimensions = H5Sget_simple_extent_ndims(space);
    if(res->num_dimensions < 0){
      sp_error_warning("Unable to get dimensions in file %s",filename);
      return NULL;
    }
    if(H5Sget_simple_extent_dims(space,dims,NULL) < 0){
      sp_error_warning("Unable to get dimensions extent in file %s",filename);
      return NULL;
    }
    if(H5Sget_simple_extent_ndims(space) == 3){
      res->image = _sp_c3matrix_alloc(dims[0],dims[1],dims[2],file,line);
      res->mask = _sp_i3matrix_alloc(dims[0],dims[1],dims[2],file,line);
    }else if(H5Sget_simple_extent_ndims(space) == 2){
      res->image = _sp_c3matrix_alloc(dims[0],dims[1],1,file,line);
      res->mask = _sp_i3matrix_alloc(dims[0],dims[1],1,file,line);
    }else{
      sp_error_warning("File has unsupported number of dimensions!\n");
      return NULL;
    }
    tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			   sp_i3matrix_z(res->mask),file,line);
    if(!tmp){
      sp_error_warning("Unable to allocate matrix");
      return NULL;
    }
    
    status = H5Dread(dataset_id,mem_type_id , H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, tmp->data);
    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }

    for(i = 0;i<sp_3matrix_size(tmp);i++){
      res->mask->data[i] = tmp->data[i];
    }
    sp_3matrix_free(tmp);
    
    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    
    dataset_id = H5Dopen(file_id, "/image_center");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }
    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);

    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }	
    

    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    res->detector->image_center[0] = values[0];
    res->detector->image_center[1] = values[1];
    if(values[2]){
      res->detector->image_center[2] = values[2];
    }else{
      res->detector->image_center[2] = 0;
    }
    
    dataset_id = H5Dopen(file_id, "/phased");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }

    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    res->phased = values[0];
    
    dataset_id = H5Dopen(file_id, "/shifted");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }

    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    res->shifted = values[0];
    
    dataset_id = H5Dopen(file_id, "/scaled");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }

    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    res->scaled = values[0];
    
    dataset_id = H5Dopen(file_id, "/detector_distance");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }

    status = H5Dread(dataset_id,  mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }

    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    res->detector->detector_distance = values[0];
    
    dataset_id = H5Dopen(file_id, "/lambda");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }

    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    res->detector->wavelength = values[0];
    
    dataset_id = H5Dopen(file_id, "/pixel_size");
    if(dataset_id < 0){
      sp_error_warning("Unable to open dataset in file %s",filename);
      return NULL;
    }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
    if(status < 0){
      sp_error_warning("Unable to read dataset from file %s",filename);
      return NULL;
    }

    status = H5Dclose(dataset_id);
    if(status < 0){
      sp_error_warning("Unable to close dataset from file %s",filename);
      return NULL;
    }

    res->detector->pixel_size[0] = values[0];
    res->detector->pixel_size[1] = values[0];
    res->detector->pixel_size[2] = values[0];

    
    if(res->phased){
      tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			     sp_i3matrix_z(res->mask),file,line);
      if(!tmp){
	sp_error_warning("Unable to allocate matrix for %s",filename);
	return NULL;
      }
      dataset_id = H5Dopen(file_id, "/complex");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }
      
      for(i = 0;i<sp_3matrix_size(tmp);i++){
	sp_imag(res->image->data[i]) = tmp->data[i];
	sp_real(res->image->data[i]) = 0;
      }
      sp_3matrix_free(tmp);
      tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			     sp_i3matrix_z(res->mask),file,line);
      if(!tmp){
	sp_error_warning("Unable to allocate matrix");
      }
      dataset_id = H5Dopen(file_id, "/real");
      if(dataset_id < 0){
	sp_error_warning("Unable to open dataset in file %s",filename);
	return NULL;
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_warning("Unable to read dataset from file %s",filename);
	return NULL;
      }

      status = H5Dclose(dataset_id);
      if(status < 0){
	sp_error_warning("Unable to close dataset from file %s",filename);
	return NULL;
      }

      for(i = 0;i<sp_3matrix_size(tmp);i++){
	sp_real(res->image->data[i]) += tmp->data[i];
      }
      sp_3matrix_free(tmp);
      
    }else{
      if(!res->scaled){
	tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			       sp_i3matrix_z(res->mask),file,line);
	if(!tmp){
	  sp_error_warning("Unable to allocate matrix");
	}
	
	dataset_id = H5Dopen(file_id, "/intensities");
	if(dataset_id < 0){
	  sp_error_warning("Unable to open dataset in file %s",filename);
	  return NULL;
	}
	
	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
	if(status < 0){
	  sp_error_warning("Unable to read dataset from file %s",filename);
	  return NULL;
	}
	
	status = H5Dclose(dataset_id);
	if(status < 0){
	  sp_error_warning("Unable to close dataset from file %s",filename);
	  return NULL;
	}
	
	for(i = 0;i<sp_3matrix_size(tmp);i++){
	  sp_real(res->image->data[i]) += tmp->data[i];
	}
	sp_3matrix_free(tmp);
      }else{
	tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
				sp_i3matrix_z(res->mask),file,line);
	if(!tmp){
	  sp_error_warning("Unable to allocate matrix");
	}

	dataset_id = H5Dopen(file_id, "/amplitudes");
	if(dataset_id < 0){
	  sp_error_warning("Unable to open dataset in file %s",filename);
	  return NULL;
	}
	
	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
	if(status < 0){
	  sp_error_warning("Unable to read dataset from file %s",filename);
	  return NULL;
	}

	status = H5Dclose(dataset_id);
	if(status < 0){
	  sp_error_warning("Unable to close dataset from file %s",filename);
	  return NULL;
	}

	for(i = 0;i<sp_3matrix_size(tmp);i++){
	  sp_real(res->image->data[i]) += tmp->data[i];
	}
	sp_3matrix_free(tmp);	 
      }
    }        
    status = H5Fclose(file_id);    
    if(status < 0){
      sp_error_warning("Unable to close file from file %s",filename);
      return NULL;
    }

  }
  H5Eset_auto(func,client_data);

  /* Due to a dataformat change when there was no change in version
     we'll have to transpose the data when there is no num_dimensions field detected.
     This is because when we changed from 2D to 3D the data changed from X changing slowest to X changing fastest
  */
  if(!flag_num_dimensions && sp_image_z(res) == 1){
    Image * tmp_img = sp_image_duplicate(res,SP_COPY_DATA|SP_COPY_MASK);
    for(int x = 0;x<sp_image_x(res);x++){
      for(int y = 0;y<sp_image_y(res);y++){
	sp_image_set(res,x,y,0,tmp_img->image->data[x*sp_image_y(res)+y]);
	sp_i3matrix_set(res->mask,x,y,0,tmp_img->mask->data[x*sp_image_y(res)+y]);
		     
      }
    }
    sp_image_free(tmp_img);
  }

  return res;
  
}

void write_cxdi(const Image * img,const char * filename){
  hsize_t  dims[3];
  hid_t dataspace_id;
  hid_t dataset_id;
  hid_t file_id;
  float * buffer = sp_malloc(sp_image_size(img)*sizeof(float));
  for(int i = 0;i<sp_image_size(img);i++){
    buffer[i] = sp_real(img->image->data[i]);
  }
  dims[0] = sp_c3matrix_x(img->image);
  dims[1] = sp_c3matrix_y(img->image);
  dims[2] = sp_c3matrix_z(img->image);

  file_id = H5Fcreate(filename,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple( 3, dims, NULL );
  H5Gcreate(file_id,"data",0);
  dataset_id = H5Dcreate(file_id, "/data/data", H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, buffer);
  H5close();
}

Image * read_cxdi(const char * filename){
  hid_t dataset_id;
  hid_t file_id;
  int status;
  file_id = H5Fopen(filename, H5F_ACC_RDONLY,H5P_DEFAULT);
  dataset_id = H5Dopen(file_id,"/data/data");
  hid_t space = H5Dget_space(dataset_id);
  if(H5Sget_simple_extent_ndims(space) == 3 ||
     H5Sget_simple_extent_ndims(space) == 2){
  }else{
    sp_error_warning("File has unsupported number of dimensions!\n");
    return NULL;
  }
  hsize_t dims[3] = {1,1,1};
  if(H5Sget_simple_extent_dims(space,dims,NULL) < 0){
    sp_error_warning("Unable to get dimensions extent in file %s",filename);
    return NULL;
  }

  Image * ret = sp_image_alloc(dims[0],dims[1],dims[2]);
  float * buffer = sp_malloc(dims[0]*dims[1]*dims[2]*sizeof(float));
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, buffer);
  for(int i = 0;i<sp_image_size(ret);i++){
    ret->image->data[i] = sp_cinit(buffer[i],0);
    ret->mask->data[i] = 1;
  }  

  H5close();
  return ret;
}


Image * read_anton_datafile(hid_t file_id,hid_t dataset_id,const char * filename){
  hid_t mem_type_id;
  int status = 0;
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }
  int nframes;
  /* read number of frames and put them together in just 1 image */
  H5Dclose(dataset_id);
  dataset_id = H5Dopen(file_id, "/data/nframes");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, &nframes);
  int total_dims[3] = {0,0,1};
  hsize_t dims[nframes][3];
  for(int i = 0;i<nframes;i++){
    dims[i][0] = 1;
    dims[i][1] = 1;
    dims[i][2] = 1;
  }
  
  real * data[nframes];
  for(int frame = 0;frame<nframes;frame++){
    char fieldname[100]; 
    sprintf(fieldname,"/data/data%i",frame);
    dataset_id = H5Dopen(file_id,fieldname);
    hid_t space = H5Dget_space(dataset_id);
    if(space < 0){
      sp_error_warning("Unable to get space in file %s",filename);
      return NULL;
    }
    if(H5Sget_simple_extent_ndims(space) == 3 ||
       H5Sget_simple_extent_ndims(space) == 2){
    }else{
      sp_error_warning("File has unsupported number of dimensions!\n");
      return NULL;
    }
    if(H5Sget_simple_extent_dims(space,dims[frame],NULL) < 0){
      sp_error_warning("Unable to get dimensions extent in file %s",filename);
      return NULL;
    }
    total_dims[0] +=dims[frame][0];
    total_dims[1] = sp_max(dims[frame][1],total_dims[1]);

    data[frame] = malloc(sizeof(real)*(dims[frame][0]*dims[frame][1]*dims[frame][2]));
#if H5_VERS_MAJOR < 2 && H5_VERS_MINOR < 8
    /* HDF5 1.6 does not support int -> double conversion so we'll need a buffer */
    if(H5Tget_class(H5Dget_type(dataset_id)) == H5T_INTEGER){
      int * buffer =  malloc(sizeof(int)*(dims[frame][0]*dims[frame][1]*dims[frame][2]));
      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, buffer);
      for(int i = 0;i<dims[frame][0]*dims[frame][1]*dims[frame][2];i++){
	data[frame][i] = buffer[i];      
      }
      free(buffer);
    }else{
      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, data[frame]);
    }
#else
    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, data[frame]);
#endif
  }
  Image * ret = sp_image_alloc(total_dims[0],total_dims[1],total_dims[2]);
  
  for(int x = 0;x<total_dims[0];x++){
    int frame = 0;
    int rel_x = x;
    while(rel_x >= dims[frame][0]){
      rel_x -= dims[frame][0];
      frame++;
    }
    for(int y = 0;y<total_dims[1];y++){
      sp_image_set(ret,x,y,0,sp_cinit(data[frame][y*dims[frame][0]+rel_x],0));
      sp_image_mask_set(ret,x,y,0,1);
    }
  }
  /* For some reason the image is upside down so we'll turn it around */
  sp_image_reflect(ret,1,SP_AXIS_X);

  return ret;  
}


Image * read_tiff(const char * filename){
  Image * out = sp_malloc(sizeof(Image));
  out->detector = sp_malloc(sizeof(Detector));
  int bpp = 4;  
  int datatype = 0;
  int width,height;
  int nstrips;
  int stripsize;
  int i;
  unsigned char * img;
  float * tmpf;
  short * tmpi;
  unsigned short * tmpui;
  unsigned char * tmpuc;

  TIFF * tif; 

  tif = TIFFOpen(filename, "r");
  if(!tif){
    return NULL;
  }
  if(TIFFGetField(tif,TIFFTAG_BITSPERSAMPLE,&bpp)){
    bpp /= 8;
  }
  if(!TIFFGetField(tif,TIFFTAG_SAMPLEFORMAT,&datatype)){
    if(bpp == 1){
      datatype = SAMPLEFORMAT_VOID;
    }else if(bpp == 2){
      datatype = SAMPLEFORMAT_UINT;
    }
  }
  
  if(!TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height)){
    perror("Could not get image height!\n");
    return NULL;
  }
  if(!TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&width)){
    perror("Could not get image width!\n");
    return NULL;
  }
  
  nstrips = TIFFNumberOfStrips(tif);
  stripsize = TIFFStripSize(tif);
  img = sp_malloc(nstrips*stripsize);
  for(i = 0;i<nstrips;i++){
    TIFFReadEncodedStrip(tif,i,img+i*stripsize,stripsize);
  }
  TIFFClose(tif);
  
  out->image = sp_c3matrix_alloc(width,height,1);
  out->mask = sp_i3matrix_alloc(width,height,1);

  if(datatype == SAMPLEFORMAT_UINT){
    tmpui = (unsigned short *)img;
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i])= tmpui[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }else if(datatype == SAMPLEFORMAT_IEEEFP){
    tmpf = (float *)img;
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i]) = tmpf[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }else if(datatype == SAMPLEFORMAT_VOID){
    tmpuc = (unsigned char *)img;
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i]) = tmpuc[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }else if(datatype == SAMPLEFORMAT_INT){
    tmpi = (short *)(img);
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i]) = tmpi[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }



  for(i = 0;i<sp_c3matrix_size(out->image);i++){
    out->mask->data[i] = 1;
  }
  sp_free(img);
  out->scaled = 0;
  out->phased = 0;
  out->shifted = 0;
  out->detector->image_center[0] = width/2;
  out->detector->image_center[1] = height/2;
  out->num_dimensions = SP_2D;
  return out;
}


void write_tiff(const Image * img,const char * filename){
  float * data;
  int nstrips;
  int stripsize;
  TIFF * tif;
  int x,y;
  int width = sp_image_x(img);
  int height = sp_image_y(img);


  tif = TIFFOpen(filename, "w");  
  nstrips = height;
  stripsize = width*sizeof(float);

  TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,width);
  TIFFSetField(tif,TIFFTAG_ROWSPERSTRIP,1);
  TIFFSetField(tif,TIFFTAG_IMAGELENGTH,height);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  data = sp_malloc(nstrips*stripsize);
  for(y = 0;y<sp_image_y(img);y++){
    for(x = 0;x<sp_image_x(img);x++){
      data[x] =sp_cabs(sp_image_get(img,x,y,0));      
    }
    TIFFWriteEncodedStrip(tif,y,data,stripsize);
  }

  /*  for(i = 0;i<nstrips;i++){
    TIFFWriteEncodedStrip(tif,i,&(data[i*stripsize/4]),stripsize);
    }*/
  TIFFClose(tif);
  sp_free(data);
}


/*! Write an image to CSV format 
 */
void write_csv(const Image * img,const char * filename){
  FILE * f = fopen(filename,"w");
  int x,y;
  if(!f){
    perror("Could not write CSV");
    exit(0);
  }
  fprintf(f,"x,y,amplitude,phase,real,imaginary\n");
  for(y = 0;y<sp_image_y(img);y++){
    for(x = 0;x<sp_image_x(img);x++){
      fprintf(f,"%d,%d,%f,%f,%f,%f\n",x,y,sp_cabs(sp_image_get(img,x,y,0)),sp_carg(sp_image_get(img,x,y,0)),sp_real(sp_image_get(img,x,y,0)),sp_imag(sp_image_get(img,x,y,0)));
    }
  }
  fclose(f);
}


int write_mask_to_png(const Image * img, char * filename, int color){
  Image  * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  int ret = 1;
  int i;
  if(sp_i3matrix_z(img->mask) != 1){
    sp_image_free(res);
    fprintf(stderr,"Can't write 3D mask to png");
  }else{
    for(i = 0;i<sp_image_size(img);i++){
      res->image->data[i] = sp_cinit(res->mask->data[i],0);
    }
    ret = write_png(res,filename,color);
    sp_image_free(res);
  }
  return ret;
}


#ifdef _WIN32
  /* png_init_io seems to crash in windows using GnuWin32 libpng-1.2.8*/

#  define READFILE(file, data, length, check) \
     check=(png_size_t)fread(data,(png_size_t)1,length,file)
#  define WRITEFILE(file, data, length, check) \
     check=(png_size_t)fwrite(data,(png_size_t)1, length, file)
#define NEAR_BUF_SIZE 1024

static void
pngtest_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
   png_uint_32 check;

   WRITEFILE((FILE *)png_ptr->io_ptr,  data, length, check);
   if (check != length)
   {
      png_error(png_ptr, "Write Error");
   }
}

static void
pngtest_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
   png_size_t check;

   /* fread() returns 0 on error, so it is OK to store this in a png_size_t
    * instead of an int, which is what fread() actually returns.
    */
   READFILE((png_FILE_p)png_ptr->io_ptr, data, length, check);

   if (check != length)
   {
      png_error(png_ptr, "Read Error!");
   }
}
#endif 

Image * read_png(const char * filename){
 FILE *fp = fopen(filename, "rb");
 int i,j;
 png_uint_32 width,height;
 int bit_depth,color_type,interlace_type,compression_type,filter_method;
 png_structp png_ptr = png_create_read_struct
   (PNG_LIBPNG_VER_STRING, (png_voidp)NULL,
    NULL,NULL);
 png_infop info_ptr = png_create_info_struct(png_ptr);
 png_byte ** row_pointers;
 Image * res;
#ifdef _WIN32
 png_set_read_fn(png_ptr, (png_voidp)fp, pngtest_read_data);
#else
 png_init_io(png_ptr, fp);
#endif
 png_read_info(png_ptr, info_ptr);
 png_get_IHDR(png_ptr, info_ptr, &width, &height,
	      &bit_depth, &color_type, &interlace_type,
	      &compression_type, &filter_method);
 
 png_textp text;
 int num_text;
 if( png_get_text( png_ptr, info_ptr, &text, &num_text ) ){
   for( i=0 ; i<num_text ; i++ ){
     //     printf( "%s: %s\n", text[i].key, text[i].text);
   }
 }

 if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
   bit_depth *= 4;
 }else if(color_type == PNG_COLOR_TYPE_GRAY_ALPHA){
   bit_depth *= 2;
 }else if(color_type == PNG_COLOR_TYPE_RGB){
   bit_depth *= 3;
 }
 row_pointers = sp_malloc(sizeof(png_byte *)*height);
 for(i = 0;i<height;i++){
   row_pointers[i] = sp_malloc(sizeof(png_byte)*width*bit_depth/8);
 }
 png_read_image(png_ptr, row_pointers);
 res = sp_image_alloc(width,height,1);
 for(i = 0;i<height;i++){
   for(j = 0;j<width;j++){
     res->image->data[i*width+j] = sp_cinit(row_pointers[i][(int)(j*bit_depth/8)],0);
     res->mask->data[i*width+j] = 1;
   }
 }
 for(i = 0;i<height;i++){
   sp_free(row_pointers[i]);
 }
 sp_free(row_pointers);
 png_destroy_read_struct(&png_ptr,&info_ptr, NULL);
 return res;
}





int write_png(const Image * img,const char * filename, int color){

  if(img->num_dimensions != SP_2D){
    fprintf(stderr,"Can only write png of 2D images in write_png!\n");
    abort();
  }

  FILE *fp = fopen(filename, "wb");
  
  png_structp png_ptr; 
  png_infop info_ptr;
  int bit_depth = 8;
  int color_type;
  int interlace_type = PNG_INTERLACE_NONE;
  int compression_type = PNG_COMPRESSION_TYPE_DEFAULT;
  int filter_method = PNG_FILTER_TYPE_DEFAULT;
  int png_transforms = PNG_TRANSFORM_IDENTITY/*|PNG_TRANSFORM_INVERT_MONO*/;
  int pixel_size = 0;
  int i,x,y;
  sp_rgb color_table[256];
  real max_v,min_v;
  png_byte ** row_pointers;

  /*fclose(fp);
  return 0;*/
  max_v = 0;
  min_v = REAL_MAX;

/* Fill color tables */
  sp_colormap_create_table(color_table,color);

  if (!fp){
    perror("Couldn't open file!\n");
    abort();
    return -1;
  }
  png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, (png_voidp)NULL/*user_error_ptr*/,
     NULL/*user_error_fn*/, NULL/*user_warning_fn*/);
  if (!png_ptr){
    perror("Couldn't allocate write structure!\n");
    abort();
    return -1;
  }
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr){
    png_destroy_write_struct(&png_ptr,
			     (png_infopp)NULL);
    perror("Couldn't allocate info structure!\n");
    abort();
    return (-1);
  }
  if (setjmp(png_jmpbuf(png_ptr))){
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    perror("Couldn't setjmp!\n");
    abort();
    return (-1);
  }
  #ifdef _WIN32
    png_set_write_fn(png_ptr, (png_voidp)fp,  pngtest_write_data,NULL);
  #else
   png_init_io(png_ptr, fp);
  #endif

  
  color_type = PNG_COLOR_TYPE_RGB;
  /* 8 bits 3 channels */
  pixel_size = 3*1;
   /* png_set_compression_level(png_ptr,Z_BEST_COMPRESSION); */
  png_set_IHDR(png_ptr, info_ptr, sp_c3matrix_x(img->image), sp_c3matrix_y(img->image),
	       bit_depth, color_type, interlace_type,
	       compression_type, filter_method);
  
  

  row_pointers = png_malloc(png_ptr,sp_c3matrix_y(img->image)*sizeof(png_byte *));
  for (i=0; i<sp_c3matrix_y(img->image); i++){
    row_pointers[i] = png_malloc(png_ptr,sp_c3matrix_x(img->image)*
				 pixel_size*sizeof(png_byte));
  }
  
  /* We're gonna scale the image so that it fits on the 8 bits */
  min_v = sp_c3matrix_min(img->image,NULL);
  max_v = sp_c3matrix_max(img->image,NULL);

  for(y = 0;y<sp_c3matrix_y(img->image);y++){
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      sp_colormap_write_rgb(&(row_pointers[y][x*3]),img,color,
			    color_table,max_v,min_v,x,y,0,0);
    }
  }
  png_set_rows(png_ptr, info_ptr, row_pointers);
  
  png_write_png(png_ptr, info_ptr, png_transforms, NULL);
  png_write_flush(png_ptr);
  /* png_write_end(png_ptr, info_ptr);*/
  for(i=0; i<sp_c3matrix_y(img->image); i++){
    png_free(png_ptr,row_pointers[i]);
  }
  png_free(png_ptr,row_pointers);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fflush(fp);
  fclose(fp);
  return 0;
}

static Image * read_smv(const char * filename){
  Image * res = NULL;
  FILE * fp = fopen(filename,"rb");
  char buffer[1024];
  int header_size = 0;
  int x_size = 0;
  int y_size = 0;
  while(fgets(buffer,1024,fp)){
    /* stop loop when we find a line with original folder */
    if(strstr(buffer,"OriginalFolder")){
      break;
    }
    char * p;
    if(strstr(buffer,"HEADER_BYTES=")){
      p = strstr(buffer,"HEADER_BYTES=")+strlen("HEADER_BYTES=");
      header_size = atoi(p);
    }
    if(strstr(buffer,"SIZE1=")){
      p = strstr(buffer,"SIZE1=")+strlen("SIZE1=");
	x_size = atoi(p);
    }
    if(strstr(buffer,"SIZE2=")){
      p = strstr(buffer,"SIZE2=")+strlen("SIZE2=");
      y_size = atoi(p);
    }
  }
  if(!x_size || !y_size || !header_size){
    return NULL;
  }
  res = sp_image_alloc(x_size,y_size,1);
  fseek(fp,header_size,SEEK_SET);
  unsigned short * data = sp_malloc(sizeof(unsigned short)*x_size*y_size);
  fread((void *)data,sizeof(unsigned short),x_size*y_size,fp);  
  for(int x = 0; x < x_size;x++){
    for(int y = 0; y < y_size;y++){
      sp_image_set(res,x,y,0,sp_cinit(data[y*x_size+x],0));
      sp_i3matrix_set(res->mask,x,y,0,1);
    }
  }
  return res;
}

/* Superseeded by the new write_vtk */
/*
int write_vtk(Image * img, const char * filename){
  FILE * f = fopen(filename,"w");
  int x,y;
  if(!f){
    perror("Bad file in write_vtk!");
    abort();
  }
  fprintf(f,"# vtk DataFile Version 2.0\n");
  fprintf(f,"Generated by image_util write_vtk()\n");
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET STRUCTURED_POINTS\n");
  fprintf(f,"DIMENSIONS %d %d 1\n",sp_c3matrix_x(img->image),sp_c3matrix_y(img->image));
  fprintf(f,"ORIGIN 0 %d 0\n",sp_c3matrix_y(img->image));
  fprintf(f,"SPACING 1 -1 1\n");
  fprintf(f,"POINT_DATA %lld\n",sp_image_size(img));
  fprintf(f,"SCALARS amplitudes float 1\n");
  fprintf(f,"LOOKUP_TABLE default\n");
  fprintf(f,"%6g",sp_cabs(img->image->data[0]));
  for(y = 0;y<sp_c3matrix_y(img->image);y++){
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      fprintf(f," %g",sp_cabs(img->image->data[y*sp_c3matrix_x(img->image)+x]));    
    }
  }
*/
/*  for(i = 1;i<sp_image_size(img);i++){
    fprintf(f," %g",img->image->data[i]);    
  }*/
/*
  fprintf(f,"\n");
  fflush(f);
  fclose(f);
  return 0;
}
*/

int write_vtk(const Image * img, const char * filename){
  FILE * f = fopen(filename,"w");
  int x,y,z;
  if(!f){
    perror("Bad file in write_vtk!");
    abort();
  }
  fprintf(f,"# vtk DataFile Version 2.0\n");
  fprintf(f,"Generated by image_util write_vtk()\n");
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET STRUCTURED_POINTS\n");
  fprintf(f,"DIMENSIONS %d %d %d\n",sp_c3matrix_x(img->image),
	  sp_c3matrix_y(img->image),sp_c3matrix_z(img->image));
  fprintf(f,"ORIGIN 0 0 0\n");//changed from 0 y 0
  fprintf(f,"SPACING 1 1 1\n");//changed from 1 -1 1 when going to 3d ??
  fprintf(f,"POINT_DATA %lld\n",sp_image_size(img));
  fprintf(f,"SCALARS amplitudes float 1\n");
  fprintf(f,"LOOKUP_TABLE default\n");
  fprintf(f,"%6g",sp_cabs(img->image->data[0]));
  for(z = 0;z<sp_c3matrix_z(img->image);z++){
    for(y = 0;y<sp_c3matrix_y(img->image);y++){
      for(x = 0;x<sp_c3matrix_x(img->image);x++){
	fprintf(f," %g",sp_cabs(img->image->data[z*sp_c3matrix_x(img->image)*sp_c3matrix_y(img->image)+y*sp_c3matrix_x(img->image)+x]));
      }
    }
  }
/*  for(i = 1;i<sp_image_size(img);i++){
    fprintf(f," %g",img->image->data[i]);    
  }*/
  fprintf(f,"\n");
  fflush(f);
  fclose(f);
  return 0;
}

int write_xplor(const Image * img, const char * filename){
  FILE * f = fopen(filename,"w");
  int x,y,z;
  if(!f){
    perror("Bad file in write_xplor!");
    abort();
  }
  fprintf(f,"       3 !NTITLE\n");
  fprintf(f,"XPLOR 3D electron density map\n");
  fprintf(f,"Generated by image_util write_xplor compiled on %s\n",__DATE__);
  time_t date = time(NULL);
  fprintf(f,"File created on: %s\n",ctime(&date));
  fprintf(f,"%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  sp_c3matrix_x(img->image),0,sp_c3matrix_x(img->image)-1,
	  sp_c3matrix_y(img->image),0,sp_c3matrix_y(img->image)-1,
	  sp_c3matrix_z(img->image),0,sp_c3matrix_z(img->image)-1);
  fprintf(f,"%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\n",(double)sp_c3matrix_x(img->image),
	  (double)sp_c3matrix_y(img->image),(double)sp_c3matrix_z(img->image),90.0,90.0,90.0);
  fprintf(f,"ZYX\n");
  real avg = 0;
  for(z = 0;z<sp_c3matrix_z(img->image);z++){
    fprintf(f,"%8d\n",z);
    int newline_counter = 0;
    for(y = 0;y<sp_c3matrix_y(img->image);y++){
      for(x = 0;x<sp_c3matrix_x(img->image);x++){
	avg += sp_cabs(sp_image_get(img,x,y,z));
	fprintf(f,"%12.5e",sp_cabs(sp_image_get(img,x,y,z)));
	newline_counter++;
	if(newline_counter == 6){
	  fprintf(f,"\n");
	  newline_counter = 0;
	}
      }
    }
    /* If necessary print last newline */
    if(newline_counter){
      fprintf(f,"\n");
    }
  }
    
  fprintf(f,"%8d\n",-9999);    
  avg /= sp_image_size(img);
  real std_dev = 0;
  for(z = 0;z<sp_c3matrix_z(img->image);z++){
    for(y = 0;y<sp_c3matrix_y(img->image);y++){
      for(x = 0;x<sp_c3matrix_x(img->image);x++){
	std_dev += (avg-sp_cabs(sp_image_get(img,x,y,z)))*
	  (avg-sp_cabs(sp_image_get(img,x,y,z)));
      }
    }
  }
  std_dev = sqrt(std_dev);
  std_dev /= sp_image_size(img);
  //  fprintf(f,"%12.4e %12.4e\n",avg,std_dev);
  fflush(f);
  fclose(f);
  return 0;
}
