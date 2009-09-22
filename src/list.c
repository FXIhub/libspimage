#include <stdlib.h>
#include <string.h>
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif
#include "spimage.h"



void sp_list_sort(sp_list * l){
  /* use bubble sort! */
  int swapped;
  do{
    swapped = 0;
    for(unsigned int i = 0;i<sp_list_size(l)-1;i++){
      if(sp_list_get(l,i) > sp_list_get(l,i+1)){
	real tmp = sp_list_get(l,i);
	sp_list_set(l,i,sp_list_get(l,i+1));
	sp_list_set(l,i+1,tmp);
	swapped = 1;
      }
    }
  }while(swapped);
}

void sp_list_union(sp_list * l, sp_list * m){
  for(unsigned int i = 0;i<sp_list_size(m);i++){
    sp_list_remove_all(l,sp_list_get(m,i));
    sp_list_append(l,sp_list_get(m,i));
  }
  sp_list_sort(l);
}

void sp_list_print(sp_list * l){
  for(unsigned int i = 0;i<sp_list_size(l);i++){
    if(i){
      printf(",");
    }
    printf("%g",sp_list_get(l,i));
  }
  printf("\n");
}
