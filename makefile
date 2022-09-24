CXX = g++
FLAGS = -O2 -std=c++14
#CFLAGS=-I

BUILD_DIR = ./build
BUILD_FRED = ./build/fred


DEPS1 = headers.hpp src/LSH/hash.hpp src/LSH/hashf.hpp src/inputHandle.hpp src/operations.hpp


OBJ_LSH = hash.o hashf.o inputHandle.o operations.o lsh.o 
OBJ_FRED = config.o curve.o frechet.o interval.o point.o simplification.o 
OBJ_CUBE = hyperc.o hashfun.o inputHandle.o operations.o hypercube.o hash.o hashf.o
OBJ_CLUSTER = initialization.o reverse.o assignment.o hash.o hashf.o inputHandle.o operations.o hyperc.o hashfun.o cluster.o update.o
OBJ_SEARCH = hash.o hashf.o inputHandle.o operations.o hyperc.o hashfun.o search.o


OUT = lsh cluster cube search

all: $(OUT)

lsh: $(OBJ_LSH) $(OBJ_FRED)
	$(CXX) $(OBJ_LSH) $(OBJ_FRED) $(FLAGS) -o lsh

cluster: $(OBJ_CLUSTER) $(OBJ_FRED)
	$(CXX) $(OBJ_CLUSTER) $(OBJ_FRED) $(FLAGS) -o cluster

cube: $(OBJ_CUBE) $(OBJ_FRED)
	$(CXX) $(OBJ_CUBE) $(OBJ_FRED) $(FLAGS) -o cube

search: $(OBJ_SEARCH) $(OBJ_FRED)
	$(CXX) $(OBJ_SEARCH) $(OBJ_FRED) $(FLAGS) -o search


##### frequently used objective files #####

inputHandle.o: src/inputHandle.cpp src/inputHandle.hpp headers.hpp src/LSH/hash.hpp
	$(CXX) $(FLAGS) -c src/inputHandle.cpp -o inputHandle.o

operations.o: src/operations.cpp src/operations.hpp headers.hpp src/LSH/hash.hpp
	$(CXX) $(FLAGS) -c src/operations.cpp -o operations.o

########################

######### fred #########
config.o: src/fred/config.cpp src/fred/config.hpp
	$(CXX) $(FLAGS) -c src/fred/config.cpp -o config.o

curve.o: src/fred/curve.cpp src/fred/curve.hpp src/fred/simplification.hpp
	$(CXX) $(FLAGS) -c src/fred/curve.cpp -o curve.o

frechet.o: src/fred/frechet.cpp src/fred/frechet.hpp
	$(CXX) $(FLAGS) -c src/fred/frechet.cpp -o frechet.o

interval.o: src/fred/interval.cpp src/fred/interval.hpp
	$(CXX) $(FLAGS) -c src/fred/interval.cpp -o interval.o

point.o: src/fred/point.cpp src/fred/point.hpp
	$(CXX) $(FLAGS) -c src/fred/point.cpp -o point.o

simplification.o: src/fred/simplification.cpp src/fred/simplification.hpp
	$(CXX) $(FLAGS) -c src/fred/simplification.cpp -o simplification.o



#########  LSH  ###############

hash.o: src/LSH/hash.cpp src/LSH/hash.hpp headers.hpp src/LSH/hashf.hpp src/operations.hpp src/inputHandle.hpp
	$(CXX) $(FLAGS) -c src/LSH/hash.cpp -o hash.o

hashf.o: src/LSH/hashf.cpp src/LSH/hashf.hpp headers.hpp src/LSH/hash.hpp src/operations.hpp
	$(CXX) $(FLAGS) -c src/LSH/hashf.cpp -o hashf.o

lsh.o: lsh.cpp headers.hpp src/LSH/hash.hpp src/inputHandle.hpp
	$(CXX) $(FLAGS) -c lsh.cpp -o lsh.o

#########  HYPERCUBE  ###############

hyperc.o: src/hypercube/hyperc.cpp src/hypercube/hyperc.hpp
	$(CXX) $(FLAGS) -c src/hypercube/hyperc.cpp -o hyperc.o

hashfun.o: src/hypercube/hashfun.cpp src/hypercube/hashfun.hpp
	$(CXX) $(FLAGS) -c src/hypercube/hashfun.cpp -o hashfun.o

hypercube.o: hypercube.cpp headers.hpp src/hypercube/hyperc.hpp src/inputHandle.hpp
	$(CXX) $(FLAGS) -c hypercube.cpp -o hypercube.o


#########  CLUSTER  ########
assignment.o: src/clustering/assignment.cpp src/clustering/assignment.hpp src/clustering/reverse.hpp src/operations.hpp src/LSH/hash.hpp
	$(CXX) $(FLAGS) -c src/clustering/assignment.cpp -o assignment.o

initialization.o: src/clustering/initialization.cpp src/clustering/initialization.hpp src/clustering/reverse.hpp src/operations.hpp src/LSH/hash.hpp
	$(CXX) $(FLAGS) -c src/clustering/initialization.cpp -o initialization.o

reverse.o: src/clustering/reverse.cpp src/clustering/reverse.hpp
	$(CXX) $(FLAGS) -c src/clustering/reverse.cpp -o reverse.o

update.o: src/clustering/update.cpp src/clustering/update.hpp src/operations.hpp src/LSH/hash.hpp src/clustering/initialization.hpp
	$(CXX) $(FLAGS) -c src/clustering/update.cpp -o update.o

cluster.o: cluster.cpp headers.hpp src/inputHandle.hpp src/LSH/hash.hpp src/clustering/initialization.hpp src/clustering/assignment.hpp src/clustering/reverse.hpp src/clustering/update.hpp
	$(CXX) $(FLAGS) -c cluster.cpp -o cluster.o


##########  SEARCH  #########
search.o: search.cpp headers.hpp src/inputHandle.hpp src/LSH/hash.hpp src/hypercube/hyperc.hpp
	$(CXX) $(FLAGS) -c search.cpp -o search.o



clean:
	rm *.o $(OUT)

clean_fred:
	rm *.o $(OUT)
