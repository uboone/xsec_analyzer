all: STV_Tools.o

%.o:%.cxx
	$(CXX) $(shell root-config --cflags --libs) -O3 -I../Utils -c -o $@ $^

clean:
	$(RM) *.o
