#ifndef _SP_MAP_H_
#define _SP_MAP_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*!
  Simple map implementation. 

  Does not allow duplicate keys.
  Insertions and lookups are O(n).
  Keys stored in ascending order.
*/
typedef struct{
  sp_list * keys;
  sp_list * values;
}sp_smap;


static inline sp_smap * sp_smap_alloc(unsigned int init_size){
  sp_smap * ret = (sp_smap *)sp_malloc(sizeof(sp_smap));
  ret->keys = sp_list_alloc(init_size);
  ret->values = sp_list_alloc(init_size);
  return ret;
}

static inline void sp_smap_free(sp_smap * m){
  sp_list_free(m->keys);
  sp_list_free(m->values);
  sp_free(m);
}

static inline void sp_smap_clear(sp_smap * m){
  sp_list_clear(m->keys);
  sp_list_clear(m->values);
}

static inline sp_list * sp_smap_get_keys(sp_smap * m){
  return m->keys;
}

static inline sp_list * sp_smap_get_values(sp_smap * m){
  return m->values;
}

static inline real sp_smap_get(sp_smap * m, real key,int * ok){
  *ok = 0;
  unsigned int i;
  for(i = 0;i<sp_list_size(m->keys);i++){
    if(sp_list_get(m->keys,i) == key){
      *ok = 1;
      return sp_list_get(m->values,i);
    }
  }
  return 0;
}

static inline void sp_smap_insert(sp_smap * m, real key,real value){
  /* calculate insertion index */
  unsigned int index;
  for(index = 0;index<sp_list_size(m->keys);index++){
    if(sp_list_get(m->keys,index) > key){
      break;
    }
  }
  sp_list_insert(m->keys,index,key);
  sp_list_insert(m->values,index,value);
}

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif

