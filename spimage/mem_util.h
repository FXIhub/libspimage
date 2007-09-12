#ifndef _MEM_UTIL_H_
#define _MEM_UTIL_H_ 1


/*! Function that makes some booking on the memory allocations
 *
 * This is an internal function. The public function is:
 * void * sp_malloc(size_t n)
 * the file and line information are then added automagically 
 * (via an evil macro)
 */
spimage_EXPORT void * _sp_malloc(size_t n, char * file, int line);
#ifdef _SP_DEBUG_MEM
#define sp_malloc(n) _sp_malloc(n,__FILE__,__LINE__)
#else
#define sp_malloc(n) malloc(n)
#endif

/*! Function that makes some booking on the memory allocations
 *
 * This is an internal function. The public function is:
 * void * sp_cmalloc(size_t nmemb,size_t size)
 * the file and line information are then added automagically 
 * (via an evil macro)
 */
spimage_EXPORT void * _sp_calloc(size_t nmemb,size_t size, char * file, int line);
#ifdef _SP_DEBUG_MEM
#define sp_calloc(nmemb,size) _sp_calloc(nmemb,size,__FILE__,__LINE__)
#else
#define sp_calloc(nmemb,size) calloc(nmemb,size)
#endif

/*! Function that makes some booking on the memory deallocations
 *
 * This is an internal function. The public function is:
 * void * sp_free(void * p)
 * the file and line information are then added automagically 
 * (via an evil macro)
 */
spimage_EXPORT void _sp_free(void *p,char *file, int line);
#ifdef _SP_DEBUG_MEM
#define sp_free(p) _sp_free(p,__FILE__,__LINE__)
#else
#define sp_free(p) free(p)
#endif

/*! Function that makes some booking on the memory allocations
 *
 * This is an internal function. The public function is:
 * void * sp_realloc(void * ptr,size_t size)
 * the file and line information are then added automagically 
 * (via an evil macro)
 */
spimage_EXPORT void * _sp_realloc(void * ptr,size_t size, char * file, int line);
#ifdef _SP_DEBUG_MEM
#define sp_realloc(ptr,size) _sp_realloc(ptr,size,__FILE__,__LINE__)
#else
#define sp_realloc(ptr,size) realloc(ptr,size)
#endif


/*! Checks any left over memory (memory leaks)
 *
 */
spimage_EXPORT void sp_malloc_finalize();

#endif 
