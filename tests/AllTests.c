#include "AllTests.h"

#include <time.h>


int RunAllTests(void)
{
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();

	//	CuSuiteAddSuite(suite, linear_alg_get_suite());
	CuSuiteAddSuite(suite, phasing_get_suite());
	//		CuSuiteAddSuite(suite, image_get_suite());
	CuSuiteAddSuite(suite, proj_get_suite());
	CuSuiteAddSuite(suite, prtf_get_suite());
	CuSuiteAddSuite(suite, container_get_suite());

	
	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
	return suite->failCount;
}

int main(void)
{
  sp_srand(time(NULL));
  
	return RunAllTests();
}
