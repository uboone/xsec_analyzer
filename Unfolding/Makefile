all: Unfolder

Unfolder: Unfolder.C 
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(XSEC_ANALYZER_DIR)/Utils/Includes $(XSEC_ANALYZER_DIR)/Bin/*.o -o $@ $^

clean:
	$(RM) Unfolder
	$(RM) *~
	$(RM) stv_root_dict_rdict.pcm
