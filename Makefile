CXX = g++
CXXFLAGS = -O3 -fopenmp

all: render

render: main.o
	$(CXX) $(CXXFLAGS) -o render main.o

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f main.o render image.png