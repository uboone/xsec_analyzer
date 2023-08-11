all: analyzer

analyzer: analyzer.C stv_root_dict.o
	make -C AnalysisEvent
	make -C Selections
	$(CXX) $(shell root-config --cflags --libs) -O3 -I. -I./Selections -I./AnalysisEvent -I./Utils ./AnalysisEvent/AnalysisEvent ./Selections/SelectionBase ./Selections/CC1mu1p0pi ./Selections/CC1mu2p0pi ./Selections/CC1muNp0pi -o $@ $^

stv_root_dict.o:
	$(RM) stv_root_dict*.*
	rootcling -f stv_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o stv_root_dict.o -c stv_root_dict.cc
	$(RM) stv_root_dict.cc

.INTERMEDIATE: stv_root_dict.o

clean:
	$(RM) analyzer stv_root_dict.o stv_root_dict_rdict.pcm
	$(RM) *~
	make clean -C AnalysisEvent
	make clean -C Selections
