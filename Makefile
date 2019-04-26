PROG = digitizer_4

all: $(PROG)

OBJ = digitizer_main_2.o digitizer_4.o
INC = -I. -I/usr/local/include -I$(OFFLINE_MAIN)/include/Event -I$(OFFLINE_MAIN)/include/
LIB = -pthread -L$(OFFLINE_MAIN)/lib -lNoRootEvent

CXX = g++
CFLAGS = -Wall -g `root-config --cflags`
LDFLAGS = -g `root-config --glibs` -Wl,-rpath,$(ROOTSYS)/lib 

$(PROG): $(OBJ)
	$(CXX) -o $(PROG) $(LDFLAGS) $(OBJ) $(LIB)

clean:
	rm -f $(PROG) *.o *~ core*

.cc.o:
	$(CXX) -c $(CFLAGS) $(INC) $<
