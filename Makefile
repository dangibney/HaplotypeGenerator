CXX = g++
CXXFLAGS = -w -ggdb3 -pthread -lm -lz -ligraph -march=native -O3 -std=c++2a
INCLUDES = -I./src/igraph-0.10.12/include -I./src/igraph-0.10.12/build/include
LIBS = -L./src/igraph-0.10.12/build/src
OBJS =  ./builds/gfa-io.o ./builds/gfa-base.o ./builds/kalloc.o

all: haplotype_generator

haplotype_generator: ./src/haplotype_generator.cpp $(OBJS)
	$(CXX) $(INCLUDES) $(LIBS)  $^ -o $@ $(CXXFLAGS)

