all: ProcessNTuples

ProcessNTuples: ProcessNTuples.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -I. -I$(STV_ANALYSIS_DIR)/Selections -I$(STV_ANALYSIS_DIR)/Utils/STV_Tools -I$(STV_ANALYSIS_DIR)/Utils/Includes $(STV_ANALYSIS_DIR)/Bin/*.o -o $@ $^

clean:
	$(RM) ProcessNTuples
	$(RM) *~
	$(RM) stv_root_dict_rdict.pcm
