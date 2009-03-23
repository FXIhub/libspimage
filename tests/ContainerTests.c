
#include "AllTests.h"

void test_sp_list_size(CuTest * tc){
  sp_list * a = sp_list_alloc(0);
  CuAssertTrue(tc,sp_list_size(a) == 0);
  sp_list_append(a,1.5);
  CuAssertTrue(tc,sp_list_size(a) == 1);
  sp_list_free(a);
}

void test_sp_list_remove(CuTest * tc){
  sp_list * a = sp_list_alloc(0);
  sp_list_append(a,1.5);
  sp_list_append(a,2);
  sp_list_append(a,3);
  sp_list_append(a,1.5);
  sp_list_append(a,2);
  CuAssertTrue(tc,sp_list_size(a) == 5);
  sp_list_remove_all(a,1.5);

  CuAssertTrue(tc,sp_list_size(a) == 3);
  CuAssertDblEquals(tc,sp_list_get(a,1),3,REAL_EPSILON);

  sp_list_remove_at(a,1);
  CuAssertDblEquals(tc,sp_list_get(a,1),2,REAL_EPSILON);
  CuAssertTrue(tc,sp_list_size(a) == 2);

  sp_list_free(a);
}

void test_sp_list_grow(CuTest * tc){
  sp_list * a = sp_list_alloc(0);
  for(int i = 0;i<10000;i++){
    sp_list_append(a,i);
  }
  CuAssertTrue(tc,sp_list_size(a) == 10000);
  sp_list_free(a);
  a = sp_list_alloc(200);
  for(int i = 0;i<10000;i++){
    sp_list_append(a,i);
  }
  CuAssertTrue(tc,sp_list_size(a) == 10000);
  sp_list_free(a);
}

void test_sp_list_sort(CuTest * tc){
  sp_list * a = sp_list_alloc(0);
  sp_list_append(a,1.5);
  sp_list_append(a,1);
  sp_list_append(a,-1);
  sp_list_sort(a);
  CuAssertDblEquals(tc,sp_list_get(a,0),-1,REAL_EPSILON);
  CuAssertDblEquals(tc,sp_list_get(a,1),1,REAL_EPSILON);
  CuAssertDblEquals(tc,sp_list_get(a,2),1.5,REAL_EPSILON);
  sp_list_free(a);
}

void test_sp_list_union(CuTest * tc){
  sp_list * a = sp_list_alloc(0);
  sp_list * b = sp_list_alloc(0);
  sp_list_append(a,1.5);
  sp_list_append(b,1);
  sp_list_append(b,2);
  sp_list_union(a,b);
  CuAssertDblEquals(tc,sp_list_get(a,0),1,REAL_EPSILON);
  CuAssertDblEquals(tc,sp_list_get(a,1),1.5,REAL_EPSILON);
  CuAssertDblEquals(tc,sp_list_get(a,2),2,REAL_EPSILON);
}


CuSuite* container_get_suite(void){
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_sp_list_size);
  SUITE_ADD_TEST(suite, test_sp_list_remove);
  SUITE_ADD_TEST(suite, test_sp_list_grow);
  SUITE_ADD_TEST(suite, test_sp_list_sort);
  SUITE_ADD_TEST(suite, test_sp_list_union);
  return suite;
}
