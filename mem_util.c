#include "spimage.h"


typedef struct sp_malloc_Node{
  struct sp_malloc_Node * previous;
  void * address;
  size_t size;
  char * allocate_from;
  struct sp_malloc_Node * next;
}Malloc_Node;

static Malloc_Node sp_malloc_head = {
  .previous = NULL,
  .address = NULL,
  .size = 0,
  .allocate_from = NULL,
  .next = NULL
};


static void sp_log(char * str);
static void sp_alloc(size_t n,char * file, int line,void * retval);
static FILE * log_file = NULL;

void sp_malloc_finalize(){
  char buffer[1024];
  int leaks = 0;
  /* check all non sp_freed pointers */
  for(Malloc_Node * p = sp_malloc_head.next;p;p = p->next){
    sprintf(buffer,"%zu bytes leaked at %p from %s",p->size,p->address,p->allocate_from);
    sp_log(buffer);
    leaks++;
  }
  sprintf(buffer,"%d leaks found!",leaks);
  sp_log(buffer);
  fclose(log_file);
}


void sp_log(char * str){
  if(!log_file){
    log_file = fopen("spimage_log.txt","w");
    if(!log_file){
      perror("Could not open log file for writing");
      abort();
    }
    atexit(sp_malloc_finalize);
  }
  fprintf(log_file,"%s\n",str);
  fflush(log_file);
}



void * _sp_realloc(void * ptr, size_t size, char * file, int line){
  char buffer[1024];
  void * retval = realloc(ptr,size);
  if (!retval){
      fprintf(stderr,"virtual memory exceeded");
      abort();
  }
  /* If ptr is NULL it's just like calling malloc */
  if(!ptr){
    sp_alloc(size,file,line,retval);
    return retval;
  }
  /* Otherwise change properties of current malloc node */
  Malloc_Node * p;
  /* find end of linked list */
  for( p = &(sp_malloc_head);p && p->address != ptr;p = p->next);
  if(!p){
    fprintf(stderr,"Could not find pointer in pointer list!");
    abort();
  }    
  sprintf(buffer,"Reallocated %p to %p freeing %zd bytes and allocating %zu bytes from %s:%d",p->address,retval,
	  p->size,size,file,line);  
  sp_log(buffer);
  p->size = size;
  p->address = retval;    
  sprintf(buffer,"%s:%d",file,line);
  free(p->allocate_from);
  p->allocate_from = malloc(sizeof(char)*(strlen(buffer)+1));
/*  strcpy(p->allocate_from,buffer);*/
  return retval;
}

void sp_alloc(size_t n, char * file, int line, void * retval){
  char buffer[1024];
  Malloc_Node * p;
  if(!log_file){
    sp_log("Initializing log system");
  }
  /* find end of linked list */
  for( p = &(sp_malloc_head);p->next;p = p->next);
  Malloc_Node * mn = malloc(sizeof(Malloc_Node));
  p->next = mn;
  mn->previous = p;
  mn->size = n;
  mn->address = retval;    
  mn->next = NULL;
  sprintf(buffer,"%zu bytes allocated at %p from %s:%d",mn->size,mn->address,file,line);  
  sp_log(buffer);
  sprintf(buffer,"%s:%d",file,line);
  mn->allocate_from = malloc(sizeof(char)*(strlen(buffer)+1));
/*  strcpy(mn->allocate_from,buffer);*/
}

void * _sp_malloc(size_t n, char * file, int line){  
  void * retval;
  if (!(retval = malloc(n))){
      fprintf(stderr,"virtual memory exceeded");
      abort();
  }
  sp_alloc(n,file,line,retval);
  return (retval);
}

void * _sp_calloc(size_t nmemb, size_t size,char * file, int line){  
  void * retval;
  if (!(retval = calloc(nmemb,size))){
      fprintf(stderr,"virtual memory exceeded");
      abort();
  }
  sp_alloc(nmemb*size,file,line,retval);
  return (retval);
}

void _sp_free(void * pointer, char * file, int line){
  Malloc_Node * p;
  char buffer[1024];
  if(!pointer){
    fprintf(stderr,"Deallocating NULL pointer!");
    abort();
  }    
  /* find pointer in allocated list */
  for(p = sp_malloc_head.next;p && (p->address != pointer);p = p->next);
  if(!p){
    fprintf(stderr,"Could not find pointer in pointer list!");
    abort();
  }
  /* head is always present so there is always previous */
  p->previous->next = p->next;
  if(p->next){    
    p->next->previous = p->previous;
  }
  sprintf(buffer,"%zu bytes sp_freed at %p from %s:%d",p->size,p->address,file,line);
  sp_log(buffer);
  free(p->allocate_from);
  free(p);  
  free(pointer);
}



