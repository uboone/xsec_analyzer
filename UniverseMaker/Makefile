all: univmake

univmake: univmake.C stv_root_dict.o
	$(CXX) $(shell root-config --cflags --libs) -O3 -I. -I./Utils -I./AnalysisEvent -o $@ $^

stv_root_dict.o:
	$(RM) stv_root_dict*.*
	rootcling -f stv_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o stv_root_dict.o -c stv_root_dict.cc
	$(RM) stv_root_dict.cc

.PHONY: clean

.INTERMEDIATE: stv_root_dict.o

clean:
	$(RM) univmake stv_root_dict.o stv_root_dict_rdict.pcm
	$(RM) *~
