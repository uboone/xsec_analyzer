all: SelectionBase CC1mu1p0pi CC1mu2p0pi CC1muNp0pi

SelectionBase: SelectionBase.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

CC1muNp0pi: CC1muNp0pi.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

CC1mu2p0pi: CC1mu2p0pi.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

CC1mu1p0pi: CC1mu1p0pi.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../ -I../Utils -I../AnalysisEvent -c -o $@ $^

clean:
	$(RM) SelectionBase CC1muNp0pi CC1mu2p0pi CC1mu1p0pi *~
