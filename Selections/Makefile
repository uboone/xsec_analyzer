#all: SelectionBase.o CC1mu1p0pi.o CC1mu2p0pi.o CC1muNp0pi.o DummySelection.o SelectionFactory.o
all: SelectionBase.o CC1muNp0pi.o DummySelection.o SelectionFactory.o
	mv *.o $(STV_ANALYSIS_DIR)/Bin

%.o:%.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(STV_ANALYSIS_DIR)/Utils/Includes -I$(STV_ANALYSIS_DIR)/Utils/STV_Tools -c -o $@ $^

clean:
	$(RM) *.o
	$(RM) *~
	$(RM) stv_root_dict_rdict.pcm
