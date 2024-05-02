all: STV_Tools.o
	mv *.o $(STV_ANALYSIS_DIR)/Bin	

%.o:%.cxx
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(STV_ANALYSIS_DIR)/Utils/Includes -c -o $@ $^

clean:
	$(RM) *.o
	$(RM) *~
