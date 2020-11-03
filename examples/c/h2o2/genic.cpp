#include <stdio.h>
#include "stdlib.h"
#include "mechanism.hpp"

#define RANDOM ((float)rand()/(float)(RAND_MAX))

int main(int argc,char* argv[])
{
  int N = atoi(argv[1]);

  double temp = 1000; // K
  double pr   = 1e5;  // Pa 

  // F order
  for(int i=0; i<N;i++) printf(" %f", RANDOM*temp);
  for(int i=0; i<N;i++) printf(" %f", pr*RANDOM);
  for(int i=0; i<N;i++)
    for(int j=0; j<NS;j++) printf(" %f", RANDOM);

  printf("\n");
}
