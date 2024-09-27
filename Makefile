# Verify that ROOT is set up before attempting a build. Don't bother with the
# check if we just want to "make clean"
ifneq ($(MAKECMDGOALS),clean)
  ROOTCONFIG := $(shell command -v root-config 2> /dev/null)
  ifndef ROOTCONFIG
    define err_message
Could not find a working ROOT installation.
Please ensure that the root-config executable is on your PATH and try again.
    endef
    $(error "$(err_message)")
  endif
endif

# Figure out the correct shared library file suffix for the current OS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHARED_LIB_SUFFIX=dylib
else ifeq ($(UNAME_S),Linux)
  SHARED_LIB_SUFFIX=so
else
  $(warning Unrecognized operating system encountered.)
  SHARED_LIB_SUFFIX=so
endif

BIN_DIR := bin
INCLUDE_DIR := include
LIB_DIR := lib

ROOT_DICTIONARY := $(LIB_DIR)/dictionaries.o
SHARED_LIB := $(LIB_DIR)/libXSecAnalyzer.$(SHARED_LIB_SUFFIX)

CXXFLAGS := $(shell root-config --cflags) -O3 -I$(INCLUDE_DIR)
LDFLAGS := $(shell root-config --libs) -L$(LIB_DIR) -lXSecAnalyzer

# Source files to use when building the main shared library
SHARED_SOURCES := $(wildcard src/binning/*.cxx)
SHARED_SOURCES += $(wildcard src/selections/*.cxx)
SHARED_SOURCES += $(wildcard src/utils/*.cxx)

# Object files that will be included in the shared library
SHARED_OBJECTS := $(SHARED_SOURCES:.cxx=.o)

# Causes GNU make to auto-delete the ROOT dictionary object file when the build
# is complete
.INTERMEDIATE: $(ROOT_DICTIONARY)

all: $(SHARED_LIB) bin/ProcessNTuples bin/univmake bin/SlicePlots \
  bin/Unfolder bin/BinScheme bin/StandaloneUnfold

$(ROOT_DICTIONARY):
	rootcling -f $(LIB_DIR)/dictionaries.cc -c LinkDef.hh
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	  -fPIC -o $@ -c $(LIB_DIR)/dictionaries.cc
	$(RM) $(LIB_DIR)/dictionaries.cc

$(SHARED_OBJECTS): %.o : %.cxx
	$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $<

$(SHARED_LIB): $(ROOT_DICTIONARY) $(SHARED_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ -fPIC -shared $^

bin/ProcessNTuples: src/app/ProcessNTuples.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/univmake: src/app/univmake.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/SlicePlots: src/app/Slice_Plots.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/Unfolder: src/app/Unfolder.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/BinScheme: src/app/binscheme.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/StandaloneUnfold: src/app/standalone_unfold.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

clean:
	$(RM) $(SHARED_LIB) $(BIN_DIR)/*
	$(RM) $(SHARED_OBJECTS)
