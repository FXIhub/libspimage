#ifndef _SPERROR_H_
#define _SPERROR_H_ 1


#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#if __STDC_VERSION__ >= 199901L
  /* Lets try our luck with variable argument macros */ 
#define  sp_error_warning(...) _sp_error_warning(__FILE__,__LINE__,__VA_ARGS__)
#define  sp_error_fatal(...) _sp_error_fatal(__FILE__,__LINE__,__VA_ARGS__)

  spimage_EXPORT void _sp_error_warning(char * file, int line, char *format, ...);
  spimage_EXPORT void _sp_error_fatal(char * file, int line, char *format, ...);

#else
  spimage_EXPORT void sp_error_warning(char *format, ...);
  spimage_EXPORT void sp_error_fatal(char *format, ...); 
#endif

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif /* _SPERROR_H_ */ 
