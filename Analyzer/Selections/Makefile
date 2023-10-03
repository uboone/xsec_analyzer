all: SelectionBase.o CC1mu1p0pi.o CC1mu2p0pi.o CC1muNp0pi.o SelectionFactory.o

SelectionBase.o: SelectionBase.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

CC1muNp0pi.o: CC1muNp0pi.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

CC1mu2p0pi.o: CC1mu2p0pi.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

CC1mu1p0pi.o: CC1mu1p0pi.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

SelectionFactory.o: SelectionFactory.h
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

clean:
	$(RM) *.o
