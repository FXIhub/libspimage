#include "AllTests.h"

#include <time.h>
#include <getopt.h>

typedef struct {
	int linear_alg;
	int phasing;
	int image;
	int proj;
	int prtf;
	int container;
	int all;
}Options;

void set_defaults(Options * opt);
Options * parse_options(int argc, char ** argv);

int RunAllTests(Options * opt)
{
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();

	if(opt->linear_alg || opt->all){
		CuSuiteAddSuite(suite, linear_alg_get_suite());
	}
	if(opt->phasing || opt->all){
		CuSuiteAddSuite(suite, phasing_get_suite());
	}
	if(opt->image || opt->all){
		CuSuiteAddSuite(suite, image_get_suite());
	}
	if(opt->proj || opt->all){
		CuSuiteAddSuite(suite, proj_get_suite());
	}
	if(opt->prtf || opt->all){
		CuSuiteAddSuite(suite, prtf_get_suite());
	}
	if(opt->container || opt->all){
		CuSuiteAddSuite(suite, container_get_suite());
	}

	
	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
	return suite->failCount;
}

Options * parse_options(int argc, char ** argv){
	int c;
  static char help_text[] = 
    "    Options description:\n\
    \n\
    -l: Run linear algebra tests\n\
    -p: Run phasing tests\n\
    -i: Run image tests\n\
    -r: Run projection tests\n\
    -t: Run PRTF tests\n\
    -c: Run container tests\n\
    -a: Run all tests\n\
    -h: print this text\n\
";
  static char optstring[] = "lpirtcah";
  Options * res = calloc(1,sizeof(Options));
  set_defaults(res);

  while(1){
    c = getopt(argc,argv,optstring);
    if(c == -1){
      break;
    }
    switch(c){
    case 'l':
      res->linear_alg = 1;
      break;
    case 'p':
      res->phasing = 1;
      break;
    case 'i':
      res->image = 1;
      break;
    case 'r':
      res->proj= 1;
      break;
    case 't':
      res->prtf = 1;
      break;
    case 'c':
      res->container = 1;
      break;
    case 'a':
      res->all = 1;
      break;
    case 'h':
      printf("%s",help_text);
      exit(0);
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  if(argc == 1){
  	printf("%s",help_text);
  	exit(0);
  }
  return res;
}

void set_defaults(Options * opt){
	opt->linear_alg = 0;
	opt->phasing = 0;
	opt->image = 0;
	opt->proj = 0;
	opt->prtf = 0;
	opt->container = 0;
	opt->all = 0;
}

int main(int argc, char ** argv)
{
	sp_srand(time(NULL));
	Options * opt = parse_options(argc, argv);
	return RunAllTests(opt);
}
