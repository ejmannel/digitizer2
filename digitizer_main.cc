#include <cstdlib>
#include <iostream>
#include <stdio.h>

int digitizer(int,  char *);

using namespace std;

int main( int argc, char* argv[] ) {

  int nevents;
  int pkt;
  
  if ( argc != 3 ) {
    cerr << "usage: digitizer pkt, filename" << endl;
    exit(1);
  } else {
    cerr << "file: " << argv[2] << endl;
  }
  
  sscanf(argv[1],"%d", &pkt);
  
  nevents = digitizer( pkt, argv[2] );
  
  return nevents;
  
}
