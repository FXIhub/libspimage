#ifndef _SP_LIST_H_
#define _SP_LIST_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


typedef struct{
  real * data;
  unsigned int used;
  unsigned int size;
}sp_list;

static inline sp_list * sp_list_alloc(unsigned int init_size){
  sp_list * ret = (sp_list *)sp_malloc(sizeof(sp_list));
  ret->data = (real *)sp_malloc(sizeof(real)*init_size);
  ret->used = 0;
  ret->size = init_size;
  return ret;
}

static inline void sp_list_free(sp_list * l){
  sp_free(l->data);
  sp_free(l);
}

static inline real sp_list_get(sp_list * l, unsigned int n){
#ifndef NDEBUG
  if(n >= l->used){
    sp_error_fatal("Trying to access list out of boundaries");
  }
#endif
  return l->data[n];
}

static inline void sp_list_set(sp_list * l,unsigned int n,real value){
#ifndef NDEBUG
  if(n >= l->used){
    sp_error_fatal("Trying to write to list out of boundaries");
  }
#endif
  l->data[n] = value;
}

static inline void sp_list_grow(sp_list * l){
  /* Growth values based on QString */
  if(l->size < 4){
    l->size = 4;
  }else if(l->size < 8){
    l->size = 8;
  }else if(l->size < 12){
    l->size = 12;
  }else if(l->size < 16){
    l->size = 16;
  }else if(l->size < 20){
    l->size = 20;
  }else if(l->size < 52){
    l->size = 52;
  }else if(l->size < 116){
    l->size = 116;
  }else if(l->size < 244){
    l->size = 244;
  }else if(l->size < 500){
    l->size = 500;
  }else if(l->size < 1012){
    l->size = 1012;
  }else if(l->size < 2036){
    l->size = 2036;
  }else if(l->size < 4084){
    l->size = 4084;
  }else if(l->size < 6132){
    l->size = 6132;
  }else if(l->size < 8180){
    l->size = 8180;
  }else if(l->size < 10228){
    l->size = 10228;
  }else if(l->size < 12276){
    l->size = 12276;
  }else if(l->size < 14324){
    l->size = 14324;
  }else if(l->size < 16372){
    l->size = 16372;
  }else{
    l->size *= 2.0;
  }
  l->data = (real *)sp_realloc(l->data,sizeof(real)*l->size);
}

static inline void sp_list_append(sp_list * l, real value){
  if(l->used == l->size){
    sp_list_grow(l);
  }
  l->data[l->used] = value;
  l->used++; 
}

static inline unsigned int sp_list_size(sp_list * l){
  return l->used;
}


static inline void sp_list_clear(sp_list * l){
  l->used = 0;
}


static inline sp_list * sp_list_duplicate(sp_list * l){
  sp_list * out = sp_list_alloc(sp_list_size(l));
  memcpy(out->data,l->data,sizeof(real)*l->used);
  out->used = l->used;
  return out;
}

static inline void sp_list_insert(sp_list * l, unsigned int n, real value){
#ifndef NDEBUG
  if(n > l->used){
    sp_error_fatal("Trying to insert in the list out of boundaries");
  }
#endif

  if(l->used == l->size){
    sp_list_grow(l);
  }
  memmove(&(l->data[n+1]),&(l->data[n]),sizeof(real)*(l->used-n));
  l->data[n] = value;
  l->used++; 
}


static inline void sp_list_remove_all(sp_list * l,  real value){
  for(int i = 0;i<l->used;i++){
    if(l->data[i] >= value-FLT_EPSILON && l->data[i] <= value+FLT_EPSILON){
      memmove(&(l->data[i]),&(l->data[i+1]),sizeof(real)*(l->used-(i+1)));
      i--;      
      l->used--;
    }
  }  
}

static inline void sp_list_remove_at(sp_list * l, int pos){
  memmove(&(l->data[pos]),&(l->data[pos+1]),sizeof(real)*(l->used-(pos+1)));
  l->used--;      
}



spimage_EXPORT void sp_list_sort(sp_list * l);
spimage_EXPORT void sp_list_union(sp_list * l,sp_list * m);
spimage_EXPORT void sp_list_print(sp_list * l);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
