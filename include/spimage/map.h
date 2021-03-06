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

  static inline sp_smap * sp_smap_create_from_pair(real key, real value){
    sp_smap * ret = sp_smap_alloc(1);
    sp_smap_insert(ret, key,value);
    return ret;
  }


static inline real sp_smap_max(sp_smap *m){
  unsigned int index;
  real value;
  real max = sp_list_get(m->values,0);
  for (index = 1; index < sp_list_size(m->values); index++) {
    value = sp_list_get(m->values,index);
    if (value > max) max = value;
  }
  return max;
}

static inline real sp_smap_min(sp_smap *m){
  unsigned int index;
  real value;
  real min = sp_list_get(m->values,0);
  for (index = 1; index < sp_list_size(m->values); index++) {
    value = sp_list_get(m->values,index);
    if (value < min) min = value;
  }
  return min;
}


  static inline real sp_smap_interpolate(sp_smap * map, real x){
  sp_list * keys = sp_smap_get_keys(map);
  sp_list * values = sp_smap_get_values(map);
  unsigned int idx = 0;
  for(idx = 0;idx<sp_list_size(keys);idx++){
    if(x < sp_list_get(keys,idx)){
      break;
    }
  }
  if(idx == 0){
    return sp_list_get(values,0);
  }
  if(idx == sp_list_size(keys)){
    return sp_list_get(values,sp_list_size(keys)-1);
  }
  /* Cubic Bezier curve taken from http://en.wikipedia.org/wiki/Bézier_curve */
  real p0y = sp_list_get(values,idx-1);
  real p0x = sp_list_get(keys,idx-1);
  real p3y = sp_list_get(values,idx);
  real p3x = sp_list_get(keys,idx);
  real t = (2 - 2*(p0x - p3x)*
	    pow(8*p0x*p3x*x - 2*p3x*pow(p0x,2) - 4*x*pow(p0x,2) + 
		2*pow(p0x,3) - 2*p0x*pow(p3x,2) - 4*x*pow(p3x,2) + 
		2*pow(p3x,3) + pow((double)pow(p0x - p3x,4)*
				   (6*p0x*p3x - 16*p0x*x - 16*p3x*x + 5*pow(p0x,2) +
				    5*pow(p3x,2) + 16*pow(x,2)),0.5),
		-0.3333333333333333) +
	    2*pow(p0x - p3x,-1)*
	    pow(8*p0x*p3x*x - 2*p3x*pow(p0x,2) - 4*x*pow(p0x,2) +
		2*pow(p0x,3) - 2*p0x*pow(p3x,2) - 4*x*pow(p3x,2) +
		2*pow(p3x,3) + pow((double)pow(p0x - p3x,4)*
				   (6*p0x*p3x - 16*p0x*x -
				    16*p3x*x + 5*pow(p0x,2) + 
				     5*pow(p3x,2) + 16*pow(x,2)),0.5)
		,0.3333333333333333)
	    )/4.;  
  return 3*p0y*t*pow(1 - t,2) + p0y*pow(1 - t,3) +
    3*p3y*(1 - t)*pow(t,2) + p3y*pow(t,3);
  }

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif

