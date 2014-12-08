/*************************************************
 * File: Sequential_KD_tree.cpp
 * Author: Kate McArdle
 *
 * Main execution file for sequential KD Tree.
 */
#include <iostream>
#include <cstdio>
#include "stdlib.h"
#include <string>
#include <sstream>
#include <vector>
#include <cstdarg>
#include <set>
#include "time.h"
#include "KDTree.h"
#include "Node.h"
#include "Point.h"
using namespace std;

// Frey: data/frey.dat
#define D 560
#define N 1965

// Mnist-test: data/mnist-test.dat
//#define D 784
//#define N 10000

/* data in file was stored in binary format. Source: Chen */
bool read_data_from_file(double **data, char *filename)
{
  FILE *fp = NULL;
  if (!(fp = fopen(filename, "rb"))) {
  	cout << filename << " didn't open" << endl;
  	return false;
  }

  int i;
  int num_in, num_total;
  for (i = 0; i < N; i++)
    {
      num_total = 0;
      while (num_total < D)
        {
          num_in = fread(data[i]+num_total, sizeof(double), D, fp);
          num_total += num_in;
        }
    }

  fclose(fp);
  return true;
}

Point<D> MakePoint(double* x) {
  Point<D> result;
  for (int i = 0; i < D; ++i) {
	  result[i] = x[i];
  }
  return result;
}

int main(int argc, char **argv) {
	char* datafile = argv[1];
	cout << "Results for " << datafile << "data:" << endl;
	int k = atoi(argv[2]);

	// Read in data:
	clock_t data_read_start = clock();
	double** data;
	data = (double**) malloc(N*sizeof(double *));
	for (int i = 0; i < N; ++i) {
		data[i] = (double*) malloc(D*sizeof(double));
	}
	if (!read_data_from_file(data, datafile)) {
		printf("error reading data!\n");
		exit(1);
	}
	clock_t data_read_end = clock();
	double data_read_time = ((double) (data_read_end - data_read_start)) / CLOCKS_PER_SEC;
	printf("time to read data: %.4f\n", data_read_time);
	// Convert to array of Points:
	clock_t points_start = clock();
	Point<D>* points;
	points = (Point<D>*) malloc(N*sizeof(Point<D>));
	for (int i = 0; i < N; ++i) {
		points[i] = MakePoint(data[i]);
	}
	clock_t points_end = clock();
	double points_time = ((double) (points_end - points_start)) / CLOCKS_PER_SEC;
	printf("time to convert to point array: %.4f\n", points_time);

	// Construct kd tree with all points at once:
	double tree_start = clock();
	KDTree<D, int> kd_all(points, N);
	double tree_end = clock();
	double tree_time = ((double) (tree_end - tree_start)) / CLOCKS_PER_SEC;
	printf("time to build kd tree with all points at once: %.4f\n", tree_time);

	// from tree kd_all:
	int** kNN_graph_all;
	kNN_graph_all = (int**) malloc(N*sizeof(int*));
	for (int i = 0; i < N; ++i) {
		kNN_graph_all[i] = (int*) malloc(k*sizeof(int));
	}
	double graph_start = clock();
	kd_all.kNN_Graph(k, kNN_graph_all);
	double graph_end = clock();
	double graph_time = ((double) (graph_end - graph_start)) / CLOCKS_PER_SEC;
	printf("time to build kNN graph from tree kd_all: %.4f\n", graph_time);

	// comparing graphs:
	printf("\nneighbors of node 1 = ");
	for (int i = 0; i < k; ++i) {
		cout << kNN_graph_all[0][i]+1 << " ";
	}
	printf("\nneighbors of node 2 = ");
	for (int i = 0; i < k; ++i) {
		cout << kNN_graph_all[1][i]+1 << " ";
	}
	cout << endl;

	free(data);
	free(points);
//	free(kNN_graph);
	free(kNN_graph_all);

}
