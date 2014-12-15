/********************************************************************
 * File: Galois_BallTree.cpp
 * Author: Nazneen Rajani
 *
 * Measuring performance method from Kate's code
 *
 */
#include <iostream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <algorithm>
#include <queue>
#include <map>
#include <cmath>
// for file reading:
#include <cstring>
#include <fstream>
#include <sstream>

#include "Galois/Galois.h"
#include "Galois/Graph/FirstGraph.h"
#include "Galois/Statistic.h"

using namespace std;

struct datapoint {
	double* pt;
	int idx;

	datapoint(double* pt, int idx) {
		this->pt = pt;
		this->idx = idx;
	}
};
//typedef struct datapoint datapoint;
struct CompareNodes {
	int dim;
	CompareNodes(int dim) : dim(dim) {}
/*
	bool operator() (const datapoint& i, const datapoint& j) {
		return (i.pt[dim] < j.pt[dim]);
	}
*/	bool operator() (const datapoint* i, const datapoint* j) {
                return (i->pt[dim] < j->pt[dim]);
        }
};
struct ballnode {
    vector<struct datapoint*> data; //list of nodes it owns
    double radius;
    datapoint *pivot;
    int dim;
    bool isLeftChild;
    bool isLeaf;
    ballnode(vector<struct datapoint*> data){
        this->data = data;
        this->radius=0;
        this->pivot= NULL;
    	dim = -1;
	isLeftChild = true;
	isLeaf=false;
	}
};

double getRadius(struct datapoint *target, vector<datapoint*>& pts, int D) {
    double radius =0.0;
    //struct datapoint *child1;
    for (int k=0; k<pts.size(); k++) {
        double dist = 0.0;
        for (int i = 0; i < D; ++i) {
            dist += (target->pt[i] - pts[k]->pt[i]) * (target->pt[i] - pts[k]->pt[i]);
        }
        if(sqrt(dist)>radius){
            radius = sqrt(dist);
            //child1 = pts[k];
        }
    }
    return radius;
}

typedef Galois::Graph::FirstGraph<struct ballnode *,void,true> BallTree;
typedef BallTree::GraphNode TreeNode;
typedef pair<TreeNode, TreeNode> TreeEdge;

struct P_buildTree {
	TreeNode bnode_parent;
	int kd_lo_idx;
	int kd_hi_idx;
	int dim;
	bool isLeftChild;
	bool isLeaf;

	P_buildTree(TreeNode bnode_parent, int kd_lo_idx, int kd_hi_idx, int dim, bool isLeftChild,bool isLeaf) {
		this->bnode_parent = bnode_parent;
		this->kd_lo_idx = kd_lo_idx;
		this->kd_hi_idx = kd_hi_idx;
		this->dim = dim;
		this->isLeftChild = isLeftChild;
		this->isLeaf = isLeaf;
	}
};

int find_split_dim(vector<datapoint*>& points, int lo_idx, int hi_idx, int n_dims) {
	int best_dim = 0;
	double max_var = 0;
	int n = hi_idx - lo_idx;

	for (int d = 0; d < n_dims; ++d) {
		double mean = 0;
		double sum_squares = 0;
		for (int i = lo_idx; i < hi_idx; ++i) {
			mean += points[i]->pt[d];
			sum_squares += (points[i]->pt[d])*(points[i]->pt[d]);
		}
		mean /= n;
		double var = (sum_squares/n) - mean*mean;
		if (var > max_var) {
			max_var = var;
			best_dim = d;
		}
	}
	return best_dim;
}

struct P_addTreeEdge {
	vector<datapoint*>& datapoints;
	BallTree& tree;
	//TreeNode bnode;
	int n_dims;

	P_addTreeEdge(vector<datapoint*>& datapoints, BallTree& tree, int n_dims) : datapoints(datapoints), tree(tree) {
		//this->bnodes = bnode;
		this->n_dims = n_dims;
	}

	void operator() (P_buildTree& curr, Galois::UserContext<P_buildTree>& wl) {
		// base cases:
		if (curr.kd_hi_idx == curr.kd_lo_idx) { return; }
		
		if ((curr.kd_hi_idx - curr.kd_lo_idx) == 1) {
			vector<datapoint*> child_points;
			//struct ballnode curr_bnode = tree.getData(curr, Galois::MethodFlag::NONE);
			for (int k=0; k<datapoints.size(); k++) {
                        	child_points.push_back(datapoints.at(k));
                	}
			struct ballnode *child = new ballnode(child_points);
			child->pivot = datapoints[0];
			child->isLeaf=true;
			child->isLeftChild=curr.isLeftChild;
			TreeNode child_bnode;
			child_bnode = tree.createNode(child);	
			tree.addNode(child_bnode, Galois::MethodFlag::NONE);
			tree.addEdge(curr.bnode_parent, child_bnode);
			return;
		}
		
		vector<datapoint*> child_points;
    		//vector<datapoint*> child2_points;
		int split_dim = curr.dim;
		for (int k=curr.kd_lo_idx; k<curr.kd_hi_idx; k++) {
            		child_points.push_back(datapoints.at(k));
    		}
		struct ballnode *child = new ballnode(child_points);
		split_dim = find_split_dim(datapoints, curr.kd_lo_idx, curr.kd_hi_idx, n_dims);
		
		sort(datapoints.begin()+curr.kd_lo_idx, datapoints.begin()+curr.kd_hi_idx, CompareNodes(split_dim));
		int med_idx = curr.kd_lo_idx + (curr.kd_hi_idx-curr.kd_lo_idx)/2;
		
		child->pivot = datapoints[med_idx];
		child->radius = getRadius(child->pivot,child->data,n_dims);
		child->isLeaf=false;
                child->isLeftChild=curr.isLeftChild;
		TreeNode child_bnode;
                child_bnode = tree.createNode(child);
                tree.addNode(child_bnode, Galois::MethodFlag::NONE);
                tree.addEdge(curr.bnode_parent, child_bnode);
		//if(med_idx-curr.kd_lo_idx > 30)
		wl.push(P_buildTree(child_bnode, curr.kd_lo_idx, med_idx, (curr.dim+1)%n_dims, true,false));
		wl.push(P_buildTree(child_bnode, med_idx+1, curr.kd_hi_idx, (curr.dim+1)%n_dims, false,false));
	}
};

class CompareDist {
public:
	bool operator() (const pair<int, double>& lhs, const pair<int, double>& rhs) const {
		if (lhs.second < rhs.second) { return true; }
		return false;
	}
};

typedef Galois::Graph::FirstGraph<struct datapoint ,double,true> KNNGraph;
typedef KNNGraph::GraphNode KNNNode;

struct P_knn {
	int k;
	int n_dims;
	KNNGraph& knnGraph;
	KNNNode& knnNodes;
	BallTree& tree;
	TreeNode& root;

	P_knn(int k, int D, BallTree& tree, TreeNode& root,KNNNode& knnNodes, KNNGraph& knnGraph) : tree(tree), root(root) ,knnNodes(knnNodes),knnGraph(knnGraph) {
		this->k = k;
		this->n_dims = D;
	}

	double getDistance(const struct datapoint key, struct datapoint* curr) {
		double dist = 0.0;
		for (int i = 0; i < n_dims; ++i) {
			dist += (key.pt[i] - curr->pt[i]) * (key.pt[i] - curr->pt[i]);
		}
		return sqrt(dist);
	}

	void search_subtree(priority_queue<pair<int, double>, vector<pair<int, double> >, CompareDist>& pq, TreeNode& curr, struct datapoint key, int level) {
		struct ballnode* curr_bnode = tree.getData(curr, Galois::MethodFlag::NONE);
		if(pq.size()==k && getDistance(key,curr_bnode->pivot)>pq.top().second)
			return;
		else if(curr_bnode->isLeaf==true){
			for (int i=0; i<curr_bnode->data.size(); i++) {
            			double dist =getDistance(key, curr_bnode->data.at(i));
            			if (pq.size()==k){
                			if(dist< pq.top().second){
                    				pq.pop();
                    				pq.emplace(make_pair(curr_bnode->data.at(i)->idx,dist));
                			}
            			}
           		 	else
                    			pq.emplace(make_pair(curr_bnode->data.at(i)->idx,dist));
        		}

		} 
		    else{
			bool isLeft_flag =false;
			double max =0.0;
			for (BallTree::edge_iterator edge : tree.out_edges(curr)) {
                                TreeNode dest = tree.getEdgeDst(edge);
				struct ballnode *child = tree.getData(dest, Galois::MethodFlag::NONE);
				double dist_1 = getDistance(key, child->pivot);
				if(dist_1>max){
					max=dist_1;
                                	if (tree.getData(dest, Galois::MethodFlag::NONE)->isLeftChild) 
						isLeft_flag=true;
					else
						isLeft_flag=false;
				}
			}
            		if(isLeft_flag==true){
				for (BallTree::edge_iterator edge : tree.out_edges(curr)) {
                                	TreeNode dest = tree.getEdgeDst(edge);
                                	if (tree.getData(dest, Galois::MethodFlag::NONE)->isLeftChild) {
						search_subtree(pq, dest, key, (level+1)%n_dims);
						break;
					}
				}
				for (BallTree::edge_iterator edge : tree.out_edges(curr)) {
                                        TreeNode dest = tree.getEdgeDst(edge);
                                        if (!tree.getData(dest, Galois::MethodFlag::NONE)->isLeftChild) {
                                                search_subtree(pq, dest, key, (level+1)%n_dims);
                                                break;
                                        }
                                }
			}
			else{
				for (BallTree::edge_iterator edge : tree.out_edges(curr)) {
                                        TreeNode dest = tree.getEdgeDst(edge);
                                        if (!tree.getData(dest, Galois::MethodFlag::NONE)->isLeftChild) {
                                                search_subtree(pq, dest, key, (level+1)%n_dims);
                                                break;
                                        }
                                }
                                for (BallTree::edge_iterator edge : tree.out_edges(curr)) {
                                        TreeNode dest = tree.getEdgeDst(edge);
                                        if (tree.getData(dest, Galois::MethodFlag::NONE)->isLeftChild) {
                                                search_subtree(pq, dest, key, (level+1)%n_dims);
                                                break;
                                        }
                                }
			}
		}
	}
	void operator() (const KNNNode& node) {
		priority_queue<pair<int, double>, vector<pair<int, double> >, CompareDist> pq;
		struct datapoint key = knnGraph.getData(node, Galois::MethodFlag::NONE);
		search_subtree(pq, root, key, 0);
	}
};

bool read_data_from_binary_double(double **data, char *filename, int N, int D) {
  FILE *fp = NULL;
  if (!(fp = fopen(filename, "rb"))) {
  	cout << filename << " didn't open" << endl;
  	return false;
  }

  int i;
  int num_in, num_total;
  for (i = 0; i < N; ++i)
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

bool read_data_from_binary_uint(double **data, char *filename, int N, int D) {
  FILE *fp = NULL;
  if (!(fp = fopen(filename, "rb"))) {
  	cout << filename << " didn't open" << endl;
  	return false;
  }

	uint8_t** data_uint;
	data_uint = (uint8_t**) malloc(N*sizeof(uint8_t *));
	for (int i = 0; i < N; ++i) {
		data_uint[i] = (uint8_t*) malloc(D*sizeof(uint8_t));
	}

  int i;
  int num_in, num_total;
  for (i = 0; i < N; ++i) {
  	num_total = 0;
  	while (num_total < D) {
  		num_in = fread(data_uint[i]+num_total, sizeof(uint8_t), D, fp);
  		num_total += num_in;
  	}
  }
  fclose(fp);

  for (i = 0; i < N; ++i) {
  	for (int j = 0; j < D; ++j) {
  		data[i][j] = (double) data[i][j];
  	}
  }
  free(data_uint);
  return true;
}

bool read_data_from_libsvm(double **data, char *filename, int N, int D) {
	for (int idx = 0; idx < N; ++idx) {
		for (int dim = 0; dim < D; ++dim) {
			data[idx][dim] = 0;
		}
	}

	ifstream f(filename);
	string line;
	int data_idx = 0;
	while(getline(f, line)) {
		istringstream iss(line);
		string label = "";
		iss >> label;
		while(iss) {
			string sub = "";
			iss >> sub;
			if(sub.length() == 0) { break; }
			string feature = sub.substr(0, sub.find(':'));
			int dim = atoi(feature.c_str()) - 1; // have to decrement to be 0-based
			string val = sub.substr(sub.find(':')+1, sub.length());
			data[data_idx][dim] = atof(val.c_str());
		}
		++data_idx;
	}
	f.close();
  return true;
}

struct Times {
	double build_tree;
	double get_knn;
	double total;

	Times() {
		build_tree = 0.0;
		get_knn = 0.0;
		total = 0.0;
	}
};

enum FileFormat { BINARY_DOUBLE, BINARY_UINT8, LIBSVM };

double** read_data(char* datafile, FileFormat fileformat, int D, int N) {
	double** data;
	data = (double**) malloc(N*sizeof(double *));
	for (int i = 0; i < N; ++i) {
		data[i] = (double*) malloc(D*sizeof(double));
	}

	bool read_result = false;
	switch(fileformat) {
	case BINARY_DOUBLE :
		read_result = read_data_from_binary_double(data, datafile, N, D);
		break;
	case BINARY_UINT8 :
		read_result = read_data_from_binary_uint(data, datafile, N, D);
		break;
	case LIBSVM :
		read_result = read_data_from_libsvm(data, datafile, N, D);
		break;
	default :
		read_result = false;
	}

	if (!read_result) {
		printf("error reading data!\n");
		exit(1);
	}
	return data;
}

void run_knn(double** data, int D, int N, int k, int n_threads, Times& times) {
	/* Read in data: */
	Galois::StatManager sm;

	Galois::setActiveThreads(1);
	Galois::StatTimer totalTime("totalTime");
	totalTime.start();

	Galois::StatTimer buildTree("buildTree");
	buildTree.start();
	BallTree balltree;
	vector<struct datapoint> points; 
	points.reserve(N);
	TreeNode bnode; 
	for (int i = 0; i < N; ++i) {
		points.emplace_back(data[i], i);
	}

	/* Build KDTree */
	int split_dim = 0;


	vector<datapoint*> pts;
    	pts.reserve(N);
    	for (int j=0; j<N; j++) {
        	pts.push_back(&points.at(j));
    	}
	split_dim = find_split_dim(pts, 0, N, D);
	sort(pts.begin(), pts.end(), CompareNodes(split_dim));
        int median = N/2;
	struct ballnode *root = new ballnode(pts);
	root->pivot = &points[median]; 
	root->radius = getRadius(root->pivot,root->data,D);	
        bnode = balltree.createNode(root);
        balltree.addNode(bnode, Galois::MethodFlag::NONE);
	vector<P_buildTree> worklist;	
	worklist.emplace_back(bnode, 0, median, 1, true,false);
	worklist.emplace_back(bnode, median+1, N, 1, false,false);
	Galois::for_each(worklist.begin(), worklist.end(), P_addTreeEdge(pts, balltree, D));
	buildTree.stop();
	times.build_tree += buildTree.get();
	
	/* Build kNN Graph */
	
	Galois::setActiveThreads(n_threads);
	Galois::StatTimer knnTime("knnTime");
	knnTime.start();
	KNNGraph knnGraph;
	KNNNode* gnodes_knn;
	gnodes_knn = (KNNNode*) malloc(N*sizeof(KNNNode));
	for (int i = 0; i < N; ++i) {
		int idx = points[i].idx;
		gnodes_knn[idx] = knnGraph.createNode(points[i]);
		knnGraph.addNode(gnodes_knn[idx], Galois::MethodFlag::NONE);
	}
	Galois::do_all(knnGraph.begin(), knnGraph.end(), P_knn(k, D, balltree, bnode,gnodes_knn[median],knnGraph));
	knnTime.stop();
	times.get_knn += knnTime.get();
	totalTime.stop();
	times.total += totalTime.get();

	//free(bnode);
	//free(gnodes_knn);
	
}

struct Data {
	char* datafile;
	FileFormat fileformat;
	int D;
	int N;
	int k;

	Data(char* f, FileFormat fileformat, int D, int N, int k) : datafile(f), fileformat(fileformat), D(D), N(N), k(k) { }
};

void measure_performance() {
	int n_reps = 1;
	vector<Data> datasets;
	int n_threads[] = { 16, 8, 4, 1 };

	Data mnist_test("../../../ScalableML/data/mnist-test.dat", BINARY_DOUBLE, 784, 10000, 8);
	datasets.push_back(mnist_test);
	Data mnist_train("../../../ScalableML/data/mnist-train.dat", BINARY_DOUBLE, 784, 60000, 8);
	datasets.push_back(mnist_train);
	Data mnist8m("/scratch/02234/kmcardle/data/mnist8m", LIBSVM, 784, 8100000, 8);
	datasets.push_back(mnist8m);
	Data tinyimgs("/scratch/02234/kmcardle/data/tiny_images.bin", BINARY_UINT8, 3072, 79302017, 8);
	datasets.push_back(tinyimgs);

	Data poker("/scratch/02234/kmcardle/data/poker.t", LIBSVM, 10, 1000000, 8);
	datasets.push_back(poker);
	Data rna("/scratch/02234/kmcardle/data/cod-rna.t", LIBSVM, 8, 271617, 8);
	datasets.push_back(rna);
	Data cadata("/scratch/02234/kmcardle/data/cadata", LIBSVM, 8, 20640, 8);
	datasets.push_back(cadata);
	Data covtype("/scratch/02234/kmcardle/data/covtype.libsvm.binary", LIBSVM, 54, 581012, 8);
	datasets.push_back(covtype);
	Data year("/scratch/02234/kmcardle/data/YearPredictionMSD", LIBSVM, 90, 463715, 8);
	datasets.push_back(year);
	Data aloi("/scratch/02234/kmcardle/data/aloi", LIBSVM, 128, 108000, 8);

	for (Data dataset : datasets) {
		cout << "---------- Results for " << dataset.datafile << " data: ----------" << endl;
		double** data = read_data(dataset.datafile, dataset.fileformat, dataset.D, dataset.N);
		for (int t : n_threads) {
			printf(" +++++ %d Threads +++++\n", t);
			Times times;
			for (int rep = 0; rep < n_reps; ++rep) {
				run_knn(data, dataset.D, dataset.N, dataset.k, t, times);
			}
			printf("\n\n +++++ %d Threads +++++\n", t);
			printf("average time to build kd tree = %f seconds\n", times.build_tree/(n_reps*1000));
			printf("average time to get knn graph = %f seconds\n", times.get_knn/(n_reps*1000));
			printf("average total time            = %f seconds\n\n\n", times.total/(n_reps*1000));
		}
		free(data);
	}

}

int main(int argc, char **argv) {
	if (argc == 6) {
		char* datafile = argv[1];
		cout << "Results for " << datafile << "data:" << endl;
		int D = atoi(argv[2]);
		int N = atoi(argv[3]);
		int k = atoi(argv[4]);
		int n_threads = atoi(argv[5]);
		double** data = read_data(datafile, BINARY_DOUBLE, D, N);
		Times times;
		run_knn(data, D, N, k, n_threads, times);
		free(data);
	}
	else {
		measure_performance();
	}
}
