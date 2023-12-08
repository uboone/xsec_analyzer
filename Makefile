all: basic NTupleProcessing UniverseMaker

basic: stv_root_dict.o
	make -C Utils
	make -C Selections

NTupleProcessing: basic
	make -C NTupleProcessing

UniverseMaker: basic
	make -C UniverseMaker

#https://root-forum.cern.ch/t/cannot-find-my-rdict-pcm-files/21901
#Seems like the stv_root_dict_rdict.pcm file is expected to be within the file where the binary is built
#For now, quick hackz is to copy that file to each subdirectory that deals with saving/reading TVectors from TTree. I hate it
stv_root_dict.o:
	$(RM) stv_root_dict*.*
	rootcling -f stv_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o stv_root_dict.o -c stv_root_dict.cc
	$(RM) stv_root_dict.cc

	cp stv_root_dict.o Bin/

	cp stv_root_dict_rdict.pcm Selections/
	cp stv_root_dict_rdict.pcm NTupleProcessing/
	cp stv_root_dict_rdict.pcm UniverseMaker/
clean:
	make clean -C Bin

	make clean -C Utils
	make clean -C Selections
	make clean -C NTupleProcessing
	make clean -C UniverseMaker

	$(RM) stv_root_dict.o stv_root_dict_rdict.pcm
