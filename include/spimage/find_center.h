#ifndef _FIND_CENTER_H_
#define _FIND_CENTER_H_ 1

#include "image.h"

spimage_EXPORT Image * sp_find_center_refine_calculate_mask(Image *img, const int search_space);
spimage_EXPORT int sp_find_center_refine(Image *img, int search_radius, int crop_size, Image *precalculated_mask);
//spimage_EXPORT int sp_find_center_refine2(Image *img, int search_radius, int crop_size);
spimage_EXPORT int sp_find_center_refine_minimal_mask(Image *img, const int search_radius, int crop_radius);
#endif
