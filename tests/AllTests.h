#ifndef _ALLTESTS_H_
#define _ALLTESTS_H_ 1

#include <assert.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "spimage.h"

#include "CuTest.h"
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif

#if defined(WIN32)
#define __func__ __FUNCTION__
#endif

#define PRINT_DONE printf("%s completed successfully\n", __func__);


#define CuAssertComplexEquals(__tc,___a,___b,__delta) do{\
    Complex __a = ___a;\
    Complex __b = ___b;\
    CuAssertDblEquals(tc,sp_real(__a),sp_real(__b),__delta);\
    CuAssertDblEquals(tc,sp_imag(__a),sp_imag(__b),__delta);\
  }while(0)

    /*    CuAssertTrue(__tc,fabs(sp_real(sp_csub(__a, __b))) < __delta && fabs(sp_imag(sp_csub(__a, __b))) < __delta); \*/

CuSuite* linear_alg_get_suite();
CuSuite* image_get_suite();
CuSuite* proj_get_suite();
CuSuite* prtf_get_suite();
CuSuite* phasing_get_suite();
CuSuite* container_get_suite();

#endif
