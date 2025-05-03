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

CXXFLAGS := $(shell root-config --cflags) -I$(INCLUDE_DIR)
LDFLAGS := $(shell root-config --libs) -L$(LIB_DIR) -lXSecAnalyzer

ifneq ($(MAKECMDGOALS),debug)
  CXXFLAGS += -O3
else
  CXXFLAGS += -O0 -g
endif

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
    bin/Unfolder bin/BinScheme bin/StandaloneUnfold bin/xsroot bin/xsnotebook \
    bin/AddFakeWeights bin/AddBeamlineGeometryWeights bin/UnfolderNuMI \
    compiledb

debug: all

# Check for the presence of compiledb on the system PATH. If it is
# available, execute it to generate a compile_commands.json file
# to be used for source code indexing by clangd.
# See https://github.com/nickdiego/compiledb for more details.
ifneq (, $(shell which compiledb))
# Environment setting here suppresses an annoying warning message
compiledb:
	@env MAKEFLAGS= MFLAGS= compiledb -n make
else
compiledb:
endif

$(ROOT_DICTIONARY):
	rootcling -f $(LIB_DIR)/dictionaries.cc -c LinkDef.hh
	$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $(LIB_DIR)/dictionaries.cc
	$(RM) $(LIB_DIR)/dictionaries.cc

$(SHARED_OBJECTS): %.o : %.cxx
	$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $<

$(SHARED_LIB): $(ROOT_DICTIONARY) $(SHARED_OBJECTS)
	$(CXX) $(CXXFLAGS) $(shell root-config --libs) -o $@ -fPIC -shared $^

bin/ProcessNTuples: src/app/ProcessNTuples.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/univmake: src/app/univmake.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/SlicePlots: src/app/Slice_Plots.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/Unfolder: src/app/Unfolder.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/UnfolderNuMI: src/app/UnfolderNuMI.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/BinScheme: src/app/binscheme.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/StandaloneUnfold: src/app/standalone_unfold.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/xsroot: $(SHARED_LIB)
	cp src/app/xsroot $(BIN_DIR)

bin/xsnotebook: $(SHARED_LIB)
	cp src/app/xsnotebook $(BIN_DIR)

bin/AddFakeWeights: src/app/NuMI/addFakeWeights.cpp $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/AddBeamlineGeometryWeights: src/app/NuMI/addBeamlineGeometryWeightsToMap.cpp $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

clean:
	$(RM) $(SHARED_LIB) $(BIN_DIR)/*
	$(RM) $(SHARED_OBJECTS)
	$(RM) compile_commands.json
