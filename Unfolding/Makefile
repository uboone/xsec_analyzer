all: Unfolding test_unfolding_refactor

test_unfolding_refactor:test_unfolding_refactor.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(STV_ANALYSIS_DIR)/Utils/Includes $(STV_ANALYSIS_DIR)/Bin/*.o -o $@ $^

Unfolding: Unfolding.C 
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(STV_ANALYSIS_DIR)/Utils/Includes $(STV_ANALYSIS_DIR)/Bin/*.o -o $@ $^

clean:
	$(RM) Unfolding
	$(RM) *~
	$(RM) stv_root_dict_rdict.pcm
