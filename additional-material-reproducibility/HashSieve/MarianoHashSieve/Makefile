SHELL = /bin/sh
RM=/bin/rm -f
#CXX=g++
#CXX=kinst-ompp-papi g++
#CXX=scorep --user --nocompiler g++ -DUSING_SCOREP
CXX=icpc
DEFS=  
PROGNAME=hashsieve
INCLUDES= 
 
LIBS= -lntl -lgmp -lgsl -lgslcblas 

DEFINES=$(INCLUDES) $(DEFS)

XFLAGS = -g -O2 -msse4.2 -fopenmp -mcmodel=large -std=c++11 -Wfatal-errors #-ltcmalloc_minimal 

SRCS = hashsieve.cpp sampler.cpp
 
OBJNAMES = $(SRCS:.cpp=.o)

X_EXE = $(PROGNAME)

all: $(X_EXE)
	
%.o: %.cpp
	$(CXX) $(XFLAGS) $(DEFINES) -o $@ -c $<

$(X_EXE): $(OBJNAMES)
	$(CXX) $(XFLAGS) $(DEFINES) -o $@ $(OBJNAMES) $(LIBS)
	
	
clean:
	$(RM) *.o $(X_EXE)

test:
	echo $(OBJNAMES) 
