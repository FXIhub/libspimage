#ifndef _TIME_UTIL_H_
#define _TIME_UTIL_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*! Returns the time since the beggining of the epoch, in microseconds */
spimage_EXPORT long long int sp_gettime();

/*! Starts a timer and returns its identifying token.

  This function starts a timer and returns a token that is used to identify the timer. 
  The returned int is a token that is used in the sp_timer_stop call.
  If the return value is < 0 an error has occurred (probably no timers were available).
  At the moment the library has a limit of _sp_max_timers timers at the same time.
  _sp_max_timers is defined in time_util.c.  

  \sa sp_timer_stop()
*/
spimage_EXPORT int sp_timer_start();

/*! Stop and releases the times. Returns the elapsed time since the beggining of the timer.
  
   The token argument is required to identify the timer.
   
   \sa sp_timer_start()
 */
spimage_EXPORT long long int sp_timer_stop(int token);

/*! Returns the elapsed time since the beggining of the timer.
  
   The token argument is required to identify the timer.
   
   \sa sp_timer_start()
 */
spimage_EXPORT long long int sp_timer_elapsed(int token);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
