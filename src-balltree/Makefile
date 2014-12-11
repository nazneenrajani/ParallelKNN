CPP = g++
CFLAGS = -O2 -fopenmp -std=c++11 
LIBS = 

test: balltree.cpp
	$(CPP) $(CFLAGS)  balltree.cpp ${LIBS} -o balltree

clean: 
	rm balltree

run:
	./balltree 20640 8 /scratch/02336/nrajani/data/cadata 32 30
