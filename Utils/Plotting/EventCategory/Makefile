all: PlotEventCategory

PlotEventCategory: PlotEventCategory.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -I. -I../../ -I../../Utils -o $@ $^

clean:
	$(RM) PlotEventCategory
	$(RM) *~
