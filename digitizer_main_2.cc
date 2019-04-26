#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include "getopt.h"

void printHelp(void);

int digitizer(int,  char *, int, int);

using namespace std;

int main( int argc, char* argv[] ) {

  int nevents = 0;
  int pkt = 21351;
  int nmodules = 3;
  int nsamples = 31;
  char filename[256] = "";
  int name = 0;
  char c;
  int verbose = 0;
  
  while (( c = getopt(argc, argv, "f:p:m:t:h")) != EOF) { 
    switch (c) 
      { 
      case 'f':
	{
	  sscanf(optarg, "%s", filename);
	  name = 1;
	}
	break;
      case 'p':
	{
	  sscanf(optarg, "%d", &pkt);
	}
	break;
      case 'm':
	{
	  sscanf(optarg, "%d", &nmodules);
	}
	break;
      case 't':
	{
	  sscanf(optarg, "%d", &nsamples);
	}
	break;
      // case 'v':
      // 	{
      // 	  verbose = 1;
      // 	}
      // 	break;
      case 'h':
	{
	  printHelp();
	  exit(-1);
	}
	break;
      }
  }
  
  if (name == 0) { 
    printHelp();
    exit(-1);
  }
  
  nevents = digitizer( pkt, filename, nsamples, nmodules);
  
  return nevents;
  
}

void printHelp(void) {

  std::cout << "Useage: digitizer2 -f filename [options]" << std::endl;
  std::cout << "        -p packet ID; default = 21351" << std::endl;
  std::cout << "        -m # ADC modules; default = 3" << std::endl;
  std::cout << "        -t # time samples; default = 31" << std::endl;
  std::cout << "        -h print this message" << std::endl;

  return;
}
