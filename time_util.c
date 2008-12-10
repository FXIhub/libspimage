#include <stdlib.h>
#ifndef _WIN32
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>
#endif

#include "spimage.h"

#if defined(_MSC_VER) || defined(__MINGW32__)
#  include <time.h>
# include <windows.h>
#include <direct.h>
#ifdef _TIMEVAL_DEFINED /* also in winsock[2].h */
#define _TIMEVAL_DEFINED
struct timeval {
    long tv_sec;
    long tv_usec;
};
#endif /* _TIMEVAL_DEFINED */
#else
#  include <sys/time.h>
#endif

#if defined(_MSC_VER) || defined(__MINGW32__)
static int gettimeofday (struct timeval *tv, void* tz) 
{ 
  union { 
    long long ns100; /*time since 1 Jan 1601 in 100ns units */ 
    FILETIME ft; 
  } now; 

  GetSystemTimeAsFileTime (&now.ft); 
  tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL); 
  tv->tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL); 
  return (0); 
}
#endif

typedef struct{
  int active;
  struct timeval tv;
}sp_timer;
  

enum { _sp_max_timers = 10 };
static sp_timer timers[_sp_max_timers] = {{0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}},
					  {0,{0,0}}};

/* Returns the number of microseconds since the begging of the epoch */
long long int sp_gettime(){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec*1e6+tv.tv_usec;
}

int sp_timer_start(){
  for(int i = 0;i<_sp_max_timers;i++){
    if(timers[i].active == 0){
      timers[i].active = 1;
      gettimeofday(&(timers[i].tv),NULL);
      /* Return timer token */
      return i;
    }
  }
  /* No timers available */
  return -1;
}

/* Returns the number of microseconds since the corresponding sp_timer_start */
long long int sp_timer_stop(int token){
  int i = token;
  struct timeval tv;
  gettimeofday(&tv,NULL);
  long long int ret = (tv.tv_sec-timers[i].tv.tv_sec)*1e6 + (tv.tv_usec-timers[i].tv.tv_usec);
  timers[i].active = 0;
  return ret;
}
