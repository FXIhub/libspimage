#ifndef _IMAGE_SPHERE_H_
#define _IMAGE_SPHERE_H_

#include "image.h"
#include "linear_alg.h"

/*
 *Structure that contains the orientation of an ewald sphere.
 */

#define SP_NEAREST_KERNEL 0
#define SP_GAUSSIAN_KERNEL 1
#define SP_LINEAR_KERNEL 2
#define SP_QSPLINE_KERNEL 3
#define SP_CSPLINE_KERNEL 4
#define SP_4SPLINE_KERNEL 5

/** @defgroup direction Direction
 *  Functions that handle direction objects.
 *  @{
 */
spimage_EXPORT Rotation * sp_rot_alloc();

spimage_EXPORT void sp_rot_free(Rotation * rot);

spimage_EXPORT Rotation * sp_rot_euler(real a1, real a2, real a3);

spimage_EXPORT void sp_rot_get_euler(Rotation * rot, real *alpha, real *beta, real *gamma);

spimage_EXPORT Rotation * sp_rot_multiply(Rotation * a, Rotation * b);

spimage_EXPORT Rotation * sp_rot_transpose(Rotation * a);

spimage_EXPORT Rotation * sp_rot_inverse();

spimage_EXPORT Rotation * sp_rot_disturb(Rotation * a, real sigma);

spimage_EXPORT real sp_rot_difference(Rotation * a, Rotation * b);

//spimage_EXPORT Rotation * sp_rotation_by_axis(real v1, real v2, real v3,real angle);

spimage_EXPORT real sp_rot_determinant(Rotation * rot);

spimage_EXPORT void sp_rot_draw(Rotation * rot);

spimage_EXPORT Rotation * sp_rot_uniform();
/*@}*/



/** @defgroup Spheres Sphere
 *  Methoods for handlinge ewald spheres.
 *  @{
 */


/*! Returns the z coordinates of a ewald sphere
 *  
 * The detector_distance and pixel_size must be set in order for
 * this function to work.
 *
 */
spimage_EXPORT sp_3matrix * sp_image_sphere_z(Image * img);
/*@}*/

/*! Collects one 2D image in one single matrix on every process.
 */
//spimage_EXPORT sp_c3matrix * sp_image_allgather(Image * img);

/*! Adds the values of the 2D image ewald_xy to the 3D image img.
 */
spimage_EXPORT void sp_image_insert_ewald(Image * img, sp_3matrix * weight, Image * slice, sp_3matrix * curvature, Rotation * rot, int kernel, real sigma, int radius);

/*! Generates a simple diffraction pattern.
 */
Image * sp_image_generate_pattern(int side);


void sp_image_get_2dpattern(Image * pattern, Image * slice, sp_3matrix * curvature, Rotation * rot, int kernel);

/*! This function takes a slice of an existing 3D diffraction pattern.
 *  Make sure that the pixel size and detector distance variables are
 *  set correctly to get the curvature right.
 */
spimage_EXPORT void sp_image_get_slice(Image * space, Image * slice, sp_3matrix * slice_z, real a, Rotation * rot);

//spimage_EXPORT void sp_image_nfft_slice(Image * space, Image * slice, sp_3matrix * curvature, Rotation * rot);
#endif
