//============================================================================
// Name        : Galois_KD_tree.cpp
// Author      : Kathryn McArdle
//============================================================================

#include <iostream>
#include <cstdlib>
#include <vector>
#include <utility>

#include "Galois/Galois.h"
#include "Galois/Graph/FirstGraph.h"
#include "Galois/Statistic.h"

using namespace std;

struct KDNode {
	double* pt;
	int idx;
	int dim;
	bool isLeftChild;

	KDNode(double* pt, int idx) {
		this->pt = pt;
		this->idx = idx;
		dim = -1;
		isLeftChild = false;
	}
};

//KDNode MakeNode(double* x, int idx) {
//  KDNode result(x, idx);
//  return result;
//}

/* data in file was stored in binary format. Source: Chen */
bool read_data_from_file(double **data, char *filename, int N, int D) {
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

typedef Galois::Graph::FirstGraph<KDNode,void,true> Graph;
typedef Graph::GraphNode GNode;
typedef pair<GNode, GNode> Edge;

int main(int argc, char **argv) {
	char* datafile = argv[1];
	cout << "Results for " << datafile << "data:" << endl;
	int D = atoi(argv[2]);
	int N = atoi(argv[3]);
	int k = atoi(argv[4]);

	// Read in data:
	clock_t data_read_start = clock();
	double** data;
	data = (double**) malloc(N*sizeof(double *));
	for (int i = 0; i < N; ++i) {
		data[i] = (double*) malloc(D*sizeof(double));
	}
	if (!read_data_from_file(data, datafile, N, D)) {
		printf("error reading data!\n");
		exit(1);
	}
	clock_t data_read_end = clock();
	double data_read_time = ((double) (data_read_end - data_read_start)) / CLOCKS_PER_SEC;
	printf("time to read data: %.4f\n", data_read_time);

	Graph kdtree;
	// Convert to vector of KDNodes and insert into graph:
	clock_t points_start = clock();
	vector<KDNode> kdnodes; // vector of kdtree nodes, NOT in indexed order (starts in order, but will be resorted)
	kdnodes.reserve(N);
	GNode* gnodes; // array of graph nodes, IN indexed order: gnodes[0] corresponds to the first data point in the file
	gnodes = (GNode*) malloc(N*sizeof(GNode));
	for (int i = 0; i < N; ++i) {
		kdnodes.emplace_back(data[i], i);
		gnodes[i] = kdtree.createNode(kdnodes[i]);
		kdtree.addNode(gnodes[i]);
	}
	clock_t points_end = clock();
	double points_time = ((double) (points_end - points_start)) / CLOCKS_PER_SEC;
	printf("time to fill graph with data points (no edges): %.4f\n", points_time);


	free(data);
	free(gnodes);
}
