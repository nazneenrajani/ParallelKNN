//============================================================================
// Name        : Galois_KD_tree.cpp
// Author      : Kathryn McArdle
//============================================================================

#include <iostream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <algorithm>
#include <queue>
#include <map>
#include <cmath>

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
		isLeftChild = true;
	}
};

struct CompareNodes {
	int dim;
	CompareNodes(int dim) : dim(dim) {}

	bool operator() (const KDNode& i, const KDNode& j) {
		return (i.pt[dim] < j.pt[dim]);
	}
};

typedef Galois::Graph::FirstGraph<KDNode,void,true> KDTree;
typedef KDTree::GraphNode TreeNode;
typedef pair<TreeNode, TreeNode> TreeEdge;

struct P_buildTree {
	int gnode_parent_idx;
	int kd_lo_idx;
	int kd_hi_idx;
	int dim;
	bool isLeftChild;

	P_buildTree(int gnode_parent_idx, int kd_lo_idx, int kd_hi_idx, int dim, bool isLeftChild) {
		this->gnode_parent_idx = gnode_parent_idx;
		this->kd_lo_idx = kd_lo_idx;
		this->kd_hi_idx = kd_hi_idx;
		this->dim = dim;
		this->isLeftChild = isLeftChild;
	}
};

struct P_addTreeEdge {
	vector<KDNode>& kdnodes;
	KDTree& tree;
	TreeNode* gnodes;
	int n_dims;

	P_addTreeEdge(vector<KDNode>& kdnodes, KDTree& tree, TreeNode* gnodes, int n_dims) : kdnodes(kdnodes), tree(tree) {
		this->gnodes = gnodes;
		this->n_dims = n_dims;
	}

	void operator() (P_buildTree& curr, Galois::UserContext<P_buildTree>& wl) {
		// base cases:
		if (curr.kd_hi_idx == curr.kd_lo_idx) { return; }
		if ((curr.kd_hi_idx - curr.kd_lo_idx) == 1) {
			int gnode_idx = kdnodes[curr.kd_lo_idx].idx;
			tree.addEdge(gnodes[curr.gnode_parent_idx], gnodes[gnode_idx]);
			KDNode& n = tree.getData(gnodes[gnode_idx]);
			n.isLeftChild = curr.isLeftChild;
			n.dim = curr.dim;
			return;
		}

		sort(kdnodes.begin()+curr.kd_lo_idx, kdnodes.begin()+curr.kd_hi_idx, CompareNodes(curr.dim));
		int med_idx = curr.kd_lo_idx + (curr.kd_hi_idx-curr.kd_lo_idx)/2;
		int gnode_idx = kdnodes[med_idx].idx;
		wl.push(P_buildTree(gnode_idx, curr.kd_lo_idx, med_idx, (curr.dim+1)%n_dims, true));
		wl.push(P_buildTree(gnode_idx, med_idx+1, curr.kd_hi_idx, (curr.dim+1)%n_dims, false));
		tree.addEdge(gnodes[curr.gnode_parent_idx], gnodes[gnode_idx]);
		KDNode& n = tree.getData(gnodes[gnode_idx]);
		n.isLeftChild = curr.isLeftChild;
		n.dim = curr.dim;
	}
};

class CompareDist {
public:
	bool operator() (const pair<TreeNode&, double>& lhs, const pair<TreeNode&, double>& rhs) const {
		if (lhs.second < rhs.second) { return true; }
		return false;
	}
};

typedef Galois::Graph::FirstGraph<KDNode,double,true> KNNGraph;
typedef KNNGraph::GraphNode KNNNode;

struct P_knn {
	int k;
	int n_dims;
	KNNGraph& knnGraph;
	KNNNode* knnNodes;
	KDTree& kdtree;
	TreeNode& kdroot;

	P_knn(int k, int D, KNNGraph& knnGraph, KNNNode* knnNodes, KDTree& kdtree, TreeNode& root) : knnGraph(knnGraph), knnNodes(knnNodes), kdtree(kdtree), kdroot(root) {
		this->k = k;
		this->n_dims = D;
	}

	double getDistance(const KNNNode& k, TreeNode& c) {
//		printf("entering getDistance...\n");
		double dist = 0.0;
		KDNode& key = knnGraph.getData(k);
		KDNode& curr = kdtree.getData(c);
		for (int i = 0; i < n_dims; ++i) {
			dist += (key.pt[i] - curr.pt[i]) * (key.pt[i] - curr.pt[i]);
		}
//		printf("leaving getDistance...\n");
		return sqrt(dist);
	}

	void search_subtree(priority_queue<pair<TreeNode&, double>, vector<pair<TreeNode&, double> >, CompareDist>& pq, TreeNode& curr, const KNNNode& key, int level) {
//		printf("level = %d\n", level);
		double dist = getDistance(key, curr);
		// if pq has less than k elems in it, push curr on:
		if (pq.size() < k) {
			pair<TreeNode&, double> p(curr, dist);
			pq.push(p);
//			pq.emplace((make_pair(curr, dist)));
		}
		// otherwise, only push on if distance of curr is less than distance of max elem, and pop max elem:
		else if (dist < pq.top().second) {
			pq.pop();
			pair<TreeNode&, double> p(curr, dist);
			pq.push(p);
//			pq.emplace((make_pair(curr, dist)));
		}
		bool left_child_searched = true;
		if (knnGraph.getData(key).pt[level] < kdtree.getData(curr).pt[level]) {
			// check to see if curr has a left child. if it does, search_subtree on that child
			for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
				TreeNode dest = kdtree.getEdgeDst(edge);
				if (kdtree.getData(dest).isLeftChild) {
					search_subtree(pq, dest, key, (level+1)%n_dims);
				}
			}
		}
		else {
			left_child_searched = false;
			// check to see if curr has a right child. if it does, search_subtree on that child
			for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
				TreeNode dest = kdtree.getEdgeDst(edge);
				if (!kdtree.getData(dest).isLeftChild) {
					search_subtree(pq, dest, key, (level+1)%n_dims);
				}
			}
		}

		// as we walk back up the tree:
		if ( (pq.size() < k) || (fabs(knnGraph.getData(key).pt[level] - kdtree.getData(curr).pt[level]) < pq.top().second) ) {
			if (left_child_searched) { // search right subtree
				for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
					TreeNode dest = kdtree.getEdgeDst(edge);
					if (!kdtree.getData(dest).isLeftChild) {
						search_subtree(pq, dest, key, (level+1)%n_dims);
					}
				}
			}
			else { // search left subtree
				for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
					TreeNode dest = kdtree.getEdgeDst(edge);
					if (kdtree.getData(dest).isLeftChild) {
						search_subtree(pq, dest, key, (level+1)%n_dims);
					}
				}
			}
		}

	}

	void operator() (const KNNNode& node) {
		priority_queue<pair<TreeNode&, double>, vector<pair<TreeNode&, double> >, CompareDist> pq;
		printf("starting to search subtree\n");
		search_subtree(pq, kdroot, node, 0);
		printf("finished searching subtree\n");
		printf("size of pq = %d\n", pq.size());
		for (int i = 0; i < k; ++i) {
			printf("i = %d\n", i);
			TreeNode& kd_dest = pq.top().first;
			printf("got kd_dest from pq\n");
			KDNode& kd_node_dest = kdtree.getData(kd_dest);
			printf("got kd node from kd_dest\n");
//			for (int j = 0; j < n_dims; ++j) {
//				cout << kd_node_dest.pt[j] << " ";
//			}
//			cout << endl;
			int dest_idx = kd_node_dest.idx;
			printf("dest_idx = %d\n", dest_idx);
//			KNNNode& knn_dest = knnNodes[dest_idx];
//			printf("accessed knnNodes[dest_idx]\n");
			double dist = pq.top().second;
			printf("accessed dist\n");
			pq.pop();
			printf("popped from pq\n");
//			knnGraph.addEdge(node, knn_dest);
//			knnGraph.getEdgeData(knnGraph.addEdge(node, knn_dest)) = dist;
//			printf("added edge with value\n");
		}
	}
};

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

int main(int argc, char **argv) {
	Galois::StatManager sm;
	char* datafile = argv[1];
	cout << "Results for " << datafile << "data:" << endl;
	int D = atoi(argv[2]);
	int N = atoi(argv[3]);
	int k = atoi(argv[4]);

	/* Read in data: */
	printf("Reading in data...\n");
	Galois::StatTimer readDataTime("readData");
	readDataTime.start();
	double** data;
	data = (double**) malloc(N*sizeof(double *));
	for (int i = 0; i < N; ++i) {
		data[i] = (double*) malloc(D*sizeof(double));
	}
	if (!read_data_from_file(data, datafile, N, D)) {
		printf("error reading data!\n");
		exit(1);
	}
	readDataTime.stop();

	Galois::setActiveThreads(1);

	Galois::StatTimer totalTime("totalTime");
	totalTime.start();

	/* Convert to vector of KDNodes and insert into graph: */
	printf("building tree...\n");
	Galois::StatTimer buildTree("buildTree");
	buildTree.start();
	KDTree kdtree;
	vector<KDNode> kdnodes; // vector of kdtree nodes, NOT in indexed order (starts in order, but will be resorted)
	kdnodes.reserve(N);
	TreeNode* gnodes; // array of graph nodes, IN indexed order: gnodes[0] corresponds to the first data point in the file
	gnodes = (TreeNode*) malloc(N*sizeof(TreeNode));
	for (int i = 0; i < N; ++i) {
		kdnodes.emplace_back(data[i], i);
		gnodes[i] = kdtree.createNode(kdnodes[i]);
		kdtree.addNode(gnodes[i]);
	}

	/* Build KDTree */
	sort(kdnodes.begin(), kdnodes.end(), CompareNodes(0));
	int median = N/2;
	int root_node_idx = kdnodes[median].idx; // corresponds to the root's index in gnodes: gnodes[root_node_idx]
	KDNode& root = kdtree.getData(gnodes[root_node_idx]);
	root.dim = 0;
	vector<P_buildTree> worklist;
	worklist.emplace_back(root_node_idx, 0, median, 1, true);
	worklist.emplace_back(root_node_idx, median+1, N, 1, false);
	Galois::for_each(worklist.begin(), worklist.end(), P_addTreeEdge(kdnodes, kdtree, gnodes, D));
	buildTree.stop();

	/* Build kNN Graph */
	printf("building knn graph...\n");
	Galois::StatTimer knnTime("knnTime");
	knnTime.start();
	KNNGraph knnGraph;
	KNNNode* gnodes_knn;
	gnodes_knn = (KNNNode*) malloc(N*sizeof(KNNNode));
	for (int i = 0; i < N; ++i) {
		int idx = kdnodes[i].idx;
		if (idx == 2501) {
			printf("idx = %d\n", idx);
		}
		gnodes_knn[idx] = knnGraph.createNode(kdnodes[i]);
		knnGraph.addNode(gnodes_knn[idx]);
	}
	printf("entering galois do all...\n");
	Galois::do_all(knnGraph.begin(), knnGraph.end(), P_knn(k, D, knnGraph, gnodes_knn, kdtree, gnodes[root_node_idx]));
	knnTime.stop();

	totalTime.stop();

	free(data);
	free(gnodes);
	free(gnodes_knn);
}

//	// test graph...
//		KDNode& node_last = kdtree.getData(gnodes[N-1]);
//		for (int i = 0; i < D; ++i) {
////			cout << node_last.pt[i] << " ";
//			if (node_last.pt[i] != data[N-1][i]) {
//				cout << node_last.pt[i] << " != " << data[N-1][i] << endl;
//			}
//		}
//		cout << endl;
//		KDNode& random = kdtree.getData(gnodes[378]);
//		for (int i = 0; i < D; ++i) {
////			cout << random.pt[i] << " ";
//			if (random.pt[i] != data[378][i]) {
//				cout << random.pt[i] << " != " << data[378][i] << endl;
//			}
//		}
//		cout << endl;

//	// test tree:
//	KDNode& root_post = kdtree.getData(gnodes[root_node_idx]);
//	printf("root node = %d\n", root_post.idx);
//	for (int i = 0; i < D; ++i) {
//			cout << root_post.pt[i] << " ";
//	}
//	cout << endl;
//
//	for (auto edge : kdtree.out_edges(gnodes[root_node_idx])) {
//		GNode dest = kdtree.getEdgeDst(edge);
//		KDNode& dest_kd = kdtree.getData(dest);
//		if (dest_kd.isLeftChild) { printf("root's left child"); }
//		else { printf("root's right child"); }
//		printf("= %d\n", dest_kd.idx);
//		for (int i = 0; i < D; ++i) {
//			cout << dest_kd.pt[i] << " ";
//		}
//		cout << endl;
//	}
