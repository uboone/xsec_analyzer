all: STV_Tools.o

%.o:%.cxx
	$(CXX) $(shell root-config --cflags --libs) -O3 -c -o $@ $^

clean:
	$(RM) *.o
