all: SelectionBase.o CC1mu1p0pi.o CC1mu2p0pi.o CC1muNp0pi.o DummySelection.o SelectionFactory.o
	mv *.o $(XSEC_ANALYZER_DIR)/Bin

%.o:%.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(XSEC_ANALYZER_DIR)/Utils/Includes -I$(XSEC_ANALYZER_DIR)/Utils/STV_Tools -c -o $@ $^

clean:
	$(RM) *.o
	$(RM) *~
	$(RM) stv_root_dict_rdict.pcm
