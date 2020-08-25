MAKE=make
TARGET=extrapanneal
CXX=g++ 
CC=g++ 

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := SimulatedAnnealingExtraP
INCDIR      := SimulatedAnnealingExtraP
BUILDDIR    := obj
TARGETDIR   := bin
RESDIR      := SimulatedAnnealingExtraP/res
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o




LIBS=-lstdc++fs

#Flags, Libraries and Includes
CXXFLAGS := -std=c++17 -O3 -g -Wall -c -fmessage-length=0 -fopenmp -MMD -MP
#DEFINES := -DUSE_NAG
LIB := -fopenmp -lm -g
LDFLAGS := -fopenmp -g
INC := -I$(INCDIR) -I/usr/local/include
INCDEP := -I$(INCDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Defauilt Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
resources: directories
	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o  $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CXXFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CXXFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

#Non-File Targets
.PHONY: all remake clean cleaner resources





















SRCDIR := SimulatedAnnealingExtraP
CPP_FILES := $(wildcard SimulatedAnnealingExtraP/*.cpp SimulatedAnnealingExtraP/alglib/*.cpp)
OBJ_FILES := $(patsubst %.cpp,%.o,$(CPP_FILES))

LANG=-std=c++17
CXXFLAGS= -O3 -g -Wall -c -fmessage-length=0 -fopenmp -MMD -MP
DEFINES= -DNUSE_NAG
LDFLAGS=-fopenmp -g
LIBS=-lstdc++fs
