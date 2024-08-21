CXX = g++
CXXFLAGS = -w -ggdb3 -pthread -lm -lz -ligraph -l:libgurobi_g++8.5.a -l:libgurobi110.so -lpthread -fopenmp -march=native -O3 -std=c++2a
INCLUDES = -I./igraph-0.10.12/include -I./igraph-0.10.12/build/include -I/opt/gurobi11.0.1_linux64/gurobi1101/linux64/include
LIBS = -L./igraph-0.10.12/build/src -L/opt/gurobi11.0.1_linux64/gurobi1101/linux64/lib
OBJS =  ./builds/gfa-io.o ./builds/gfa-base.o ./builds/kalloc.o

all: ILP_constructor

ILP_constructor: ./src/ILP_constructor.cpp $(OBJS) 
	$(CXX) $(INCLUDES) $(LIBS)  $^ -o $@ $(CXXFLAGS)

ed_calc: ./src/edit_distance_calc.cpp ./src/edlib.cpp $(OBJS) 
	$(CXX) $(INCLUDES) -I./src $(LIBS)  $^ -o $@ $(CXXFLAGS)

%.o: %.cpp $(OBJS)
	$(CXX) $(INCLUDES) $(LIBS) -c $< -o $@ $(CXXFLAGS) 

haplotype_generator: ./src/haplotype_generator.cpp $(OBJS)
	$(CXX) $(INCLUDES) $(LIBS)  $^ -o $@ $(CXXFLAGS)

gurobi_solver: ./src/gurobi_solver.cpp $(OBJS)
	$(CXX) $(INCLUDES) $(LIBS)  $^ -o $@ $(CXXFLAGS)

solution_extractor: ./src/solution_extractor.cpp $(OBJS)
	$(CXX) $(INCLUDES) $(LIBS)  $^ -o $@ $(CXXFLAGS)

read_aligner: ./src/read_aligner.cpp $(OBJS) 
	$(CXX) $(INCLUDES) $(LIBS)  $^ -o $@ $(CXXFLAGS)

remove_reads_in_all: ./src/remove_reads_in_all.cpp $(OBJS) 
	$(CXX) $(INCLUDES) $(LIBS)  $^ -o $@ $(CXXFLAGS)

clean:
	rm -f *.o test_tool
