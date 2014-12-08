CC = /opt/local/bin/g++
OPENMP = -fopenmp
CFLAGS = -g -Wall -std=c++11
LIBS = 

TARGET = kdtree

all: $(TARGET)

kdtree: Sequential_KD_tree.o
	$(CC) -o $@ $(LIBS) $(OPENMP) Sequential_KD_tree.o

Sequential_KD_tree.o: Sequential_KD_tree.cpp KDTree.h Node.h Point.h
	$(CC) -c $(OPENMP) $(CFLAGS) Sequential_KD_tree.cpp

clean:
	rm -f *.o $(TARGET)

frey:
	./$(TARGET) "../data/frey.dat" 12

mnist-test:
	./$(TARGET) "../data/mnist-test.dat" 8
