all: ProcessNTuples

ProcessNTuples: ProcessNTuples.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -I. -I$(XSEC_ANALYZER_DIR)/Selections -I$(XSEC_ANALYZER_DIR)/Utils/STV_Tools -I$(XSEC_ANALYZER_DIR)/Utils/Includes $(XSEC_ANALYZER_DIR)/Bin/*.o -o $@ $^

clean:
	$(RM) ProcessNTuples
	$(RM) *~
	$(RM) stv_root_dict_rdict.pcm
