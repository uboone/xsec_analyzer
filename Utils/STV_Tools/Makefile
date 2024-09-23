all: STV_Tools.o
	mv *.o $(XSEC_ANALYZER_DIR)/Bin	

%.o:%.cxx
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(XSEC_ANALYZER_DIR)/Utils/Includes -c -o $@ $^

clean:
	$(RM) *.o
	$(RM) *~
