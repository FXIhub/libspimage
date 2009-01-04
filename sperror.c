#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

const char      lib_name[] = "libspimage";

/* Issue error message to standard error stream and, if status value not
 * less than zero, terminate calling program.  Argument mode is intended
 * to describe error severity. */

  /* Lets try our luck with variable argument macros */ 
void sp_error_report(char * file, int line, int status, char *mode, char *format,va_list ap){
   fprintf(stderr, "%s: %s: ", lib_name, mode);
   vfprintf(stderr, format, ap);
   fprintf(stderr, " in %s:%d\n",file,line);
   if (status >= 0){
     exit(status);
   }
}

void _sp_error_warning(char * file, int line, char *format, ...){
  va_list ap;
  va_start(ap,format);
  sp_error_report(file,line,-1, "warning", format,ap);
  va_end(ap);
}

void _sp_error_fatal(char * file, int line, char *format, ...){
  va_list ap;
  va_start(ap,format);
  sp_error_report(file, line,EXIT_FAILURE, "FATAL", format,ap);
  va_end(ap);
}

 void sp_error_report2(int status,  char *mode,  char *format,va_list ap){
   fprintf(stderr, "%s: %s: ", lib_name, mode);
   vfprintf(stderr, format, ap);
   fprintf(stderr, "\n");
   if (status >= 0){
     exit(status);
   }
}

void sp_error_warning( char *format, ...){
  va_list ap;
  va_start(ap,format);
  sp_error_report2(-1, "warning", format,ap);
  va_end(ap);
}

void sp_error_fatal( char *format, ...){
  va_list ap;
  va_start(ap,format);
  sp_error_report2(EXIT_FAILURE, "FATAL", format,ap);
  va_end(ap);
}

