all: SelectionBase.o CC1mu1p0pi.o CC1mu2p0pi.o CC1muNp0pi.o DummySelection.o SelectionFactory.o

%.o:%.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

%.o:%.h
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

clean:
	$(RM) *.o
