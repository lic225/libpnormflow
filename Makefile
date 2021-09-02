CXX = g++-11
CXXFLAGS = -O3 -std=c++20 -Wall -Wshadow

all: main

main: main.cpp linalg.o graph.o num_methods.o
	$(CXX) $(CXXFLAGS) -o main $^

num_methods.o: num_methods.cpp num_methods.h
	$(CXX) $(CXXFLAGS) -c -o $@ num_methods.cpp

linalg.o: linalg.cpp linalg.h
	$(CXX) $(CXXFLAGS) -c -o $@ linalg.cpp

graph.o: graph.cpp graph.h
	$(CXX) $(CXXFLAGS) -c -o $@ graph.cpp

clean:
	rm main *.o
