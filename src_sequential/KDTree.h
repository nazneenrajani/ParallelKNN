/********************************************************************
 * File: KDTree.h
 * Author: Kate McArdle
 *
 * Based on Stanford KDTree assignment: http://web.stanford.edu/class/cs106l/handouts/assignment-3-kdtree.pdf
 *
 * An interface representing a kd-tree in some number of dimensions.
 * The tree can be constructed from a set of data and then queried
 * for membership and nearest neighbors.
 */

#ifndef KDTree_Included
#define KDTree_Included

#include "Point.h"
//#include "BoundedPQueue.h"
#include "Node.h"
#include <stdexcept>
#include <cmath>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <string>
using namespace std;

template <size_t N, typename ElemType> class KDTree {
public:
	class CompareDist {
	public:
		bool operator() (const pair<int, double>& lhs, const pair<int, double>& rhs) const {
			if (lhs.second < rhs.second) { return true; }
			return false;
		}
	};
	/* Constructor: KDTree();
	 * Usage: KDTree myTree;
	 * ----------------------------------------------------
	 * Constructs an empty KDTree.
	 */
	KDTree(); // done
	void build_subtree(Node<N, ElemType>* parent, string child, vector<Node<N, ElemType>* >& nodes, int start_idx, int end_idx, int dim);
	KDTree(Point<N>* points, int n_items);

	/* Destructor: ~KDTree()
	 * Usage: (implicit)
	 * ----------------------------------------------------
	 * Cleans up all resources used by the KDTree.
	 */
	void remove_node(Node<N, ElemType>* node);
	~KDTree(); // done

	/* KDTree(const KDTree& rhs);
	 * KDTree& operator= (const KDTree& rhs);
	 * Usage: KDTree one = two;
	 * Usage: one = two;
	 * -----------------------------------------------------
	 * Deep-copies the contents of another KDTree into this
	 * one.
	 */
	KDTree(const KDTree& other); // done
	KDTree& operator= (const KDTree& other);

	/* size_t dimension() const;
	 * Usage: size_t dim = kd.dimension();
	 * ----------------------------------------------------
	 * Returns the dimension of the points stored in this
	 * KDTree.
	 */
	size_t dimension() const; // done

	/* size_t size()  const;
	 * bool   empty() const;
	 * Usage: if (kd.empty())
	 * ----------------------------------------------------
	 * Returns the number of elements in the kd-tree and
	 * whether the tree is empty.
	 */
	size_t size() const; // done
	bool   empty() const; // done

	/* bool contains(const Point<N>& pt) const;
	 * Usage: if (kd.contains(pt)) { ... }
	 * ----------------------------------------------------
	 * Returns whether the specified point is contained in 
	 * the KDTree.
	 */
	bool contains(const Point<N>& pt) const; // done

	/* void insert(const Point<N>& pt, const ElemType& value);
	 * Usage: kd.insert(v, "This value is associated with v.");
	 * ----------------------------------------------------
	 * Inserts the point pt into the KDTree, associating it
	 * with the specified value.  If the element already existed
	 * in the tree, the new value will overwrite the existing
	 * one.
	 */
	void insert(const Point<N>& pt, const ElemType& value); // done

	/* ElemType& operator[] (const Point<N>& pt);
	 * Usage: kd[v] = "Some Value";
	 * ----------------------------------------------------
	 * Returns a reference to the value associated with point
	 * pt in the KDTree.  If the point does not exist, then
	 * it is added to the KDTree using the default value of
	 * ElemType as its key.
	 */
	ElemType& operator[] (const Point<N>& pt); // done

	/* ElemType& at(const Point<N>& pt);
	 * const ElemType& at(const Point<N>& pt) const;
	 * Usage: cout << kd.at(v) << endl;
	 * ----------------------------------------------------
	 * Returns a reference to the key associated with the point
	 * pt.  If the point is not in the tree, this function throws
	 * an out_of_range exception.
	 */
	ElemType& at(const Point<N>& pt); // done
	const ElemType& at(const Point<N>& pt) const; // done

	/* ElemType kNNValue(const Point<N>& key, size_t k) const
	 * Usage: cout << kd.kNNValue(v, 3) << endl;
	 * ----------------------------------------------------
	 * Given a point v and an integer k, finds the k points
	 * in the KDTree nearest to v and returns the most common
	 * value associated with those points.  In the event of
	 * a tie, one of the most frequent value will be chosen.
	 */
	void search_subtree(priority_queue<pair<int, double>, vector<pair<int, double> >, CompareDist>& pq, const Node<N, ElemType>* curr, const Point<N>& key, int level, int k) const; // done
	ElemType kNNValue(const Point<N>& key, int k) const; // done
	vector<const Node<N, ElemType>*> get_kNN(const Point<N>& key, int k) const;
	void kNN_Graph(int k, ElemType** neighbors) const;
	void kNN_Graph_recursive(int k, ElemType** neighbors, Node<N, ElemType>* node) const;

private:
  Node<N, ElemType>* root;
  size_t n_items;
  Node<N, ElemType>* findNode(const Point<N>& pt) const; // returns a reference to the point if found; else returns NULL
};

/* * * * * Implementation Below This Point. * * * * */

template <size_t N, typename ElemType> KDTree<N, ElemType>::KDTree() {
  n_items = 0;
  root = NULL;
}

template<size_t N, typename ElemType>
struct CompareNodes {
	int dim;
	CompareNodes(int dim) { this->dim = dim; }

	bool operator() (Node<N, ElemType>* i, Node<N, ElemType>* j) {
		return (i->getPoint()[dim] < j->getPoint()[dim]);
	}
};

template<size_t N, typename ElemType> void KDTree<N, ElemType>::build_subtree(Node<N, ElemType>* parent, string child, vector<Node<N, ElemType>* >& nodes, int start_idx, int end_idx, int dim) {
	// base cases:
	if (end_idx == start_idx) { return; }

	if ((end_idx - start_idx) == 1) {
		if (child == "left") {
			parent->setLeftChild(nodes[start_idx]);
		}
		else if (child == "right") {
			parent->setRightChild(nodes[start_idx]);
		}
		return;
	}

	sort(nodes.begin()+start_idx, nodes.begin()+end_idx, CompareNodes<N, ElemType>(dim));
	int med_idx = start_idx + (end_idx-start_idx)/2;
	if (child == "left") {
		parent->setLeftChild(nodes[med_idx]);
		build_subtree(parent->getLeftChild(), "left", nodes, start_idx, med_idx, (dim+1)%N);
		build_subtree(parent->getLeftChild(), "right", nodes, med_idx+1, end_idx, (dim+1)%N);
	}
	else if (child == "right") {
		parent->setRightChild(nodes[med_idx]);
		build_subtree(parent->getRightChild(), "left", nodes, start_idx, med_idx, (dim+1)%N);
		build_subtree(parent->getRightChild(), "right", nodes, med_idx+1, end_idx, (dim+1)%N);
	}
}

template <size_t N, typename ElemType> KDTree<N, ElemType>::KDTree(Point<N>* points, int n_items) {
  this->n_items = n_items;
  // Make vector of nodes out of points
  vector<Node<N, ElemType>* > nodes;
  nodes.reserve(n_items);
  for (int i = 0; i < n_items; ++i) {
  	nodes.push_back(new Node<N, ElemType>(points[i], i));
  }

  // sort points on first dimension, get median and assign it as root
  sort(nodes.begin(), nodes.end(), CompareNodes<N, ElemType>(0));
  int med_idx = n_items/2;
  root = nodes[med_idx];

  // recursively add left and right children until all points are in tree
  build_subtree(root, "left", nodes, 0, med_idx, 1);
  build_subtree(root, "right", nodes, med_idx+1, n_items, 1);
}

template <size_t N, typename ElemType> KDTree<N, ElemType>::KDTree(const KDTree& other) {
  n_items = other.n_items;
  if (other.root != NULL) {
	  root = new Node<N, ElemType>(*(other.root));
  }
  else { root = NULL; }
}

template <size_t N, typename ElemType> KDTree<N, ElemType>& KDTree<N, ElemType>::operator= (const KDTree& other) {
  if (this != &other) {
	  n_items = other.n_items;
	  if (other.root != NULL) {
		  if (root != NULL) {
			  remove_node(root);
		  }
		  root = new Node<N, ElemType>(*(other.root));
	  }
	  else { root = NULL; }
  }
  return *this;
}

template <size_t N, typename ElemType> void KDTree<N, ElemType>::remove_node(Node<N, ElemType>* node) {
	Node<N, ElemType>* lc = node->getLeftChild();
	if (lc != NULL) remove_node(lc);
	Node<N, ElemType>* rc = node->getRightChild();
	if (rc != NULL) remove_node(rc);
	delete node;
}

template <size_t N, typename ElemType> KDTree<N, ElemType>::~KDTree() {
  if (root != NULL) {
	  remove_node(root);
  }
}



template <size_t N, typename ElemType> size_t KDTree<N, ElemType>::dimension() const {
  return N;
}

template <size_t N, typename ElemType> size_t KDTree<N, ElemType>::size() const {
	return n_items;
}

template <size_t N, typename ElemType> bool KDTree<N, ElemType>::empty() const {
	return (n_items == 0);
}

template <size_t N, typename ElemType> void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
	if (root == NULL) {
		root = new Node<N, ElemType>(pt, value);
		++n_items;
		return;
	}
	else {
		bool newNodeAdded = root->add(pt, value, 0);
		if (newNodeAdded) {
			++n_items;
		}
		return;
	}
}

template <size_t N, typename ElemType> Node<N, ElemType>* KDTree<N, ElemType>::findNode(const Point<N>& pt) const {
	Node<N, ElemType>* current = root;
	int level = 0;
	while (current != NULL) {
		if(current->getPoint() == pt) {
			return current;
		}
		else if (pt[level%N] < current->getPoint()[level%N]) {
			current = current->getLeftChild();
		}
		else {
			current = current->getRightChild();
		}
		++level;
	}
	return NULL;
}

template <size_t N, typename ElemType> bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
	return (findNode(pt) != NULL);
}

template <size_t N, typename ElemType> ElemType& KDTree<N, ElemType>::operator[] (const Point<N>& pt) {
	Node<N, ElemType>* n = findNode(pt);
	if (n != NULL) {
		return n->getVal();
	}
	else {
		ElemType* d = new ElemType();
		this->insert(pt, *d);
		n = findNode(pt);
		return n->getVal();
	}
}

template <size_t N, typename ElemType> ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
	Node<N, ElemType>* n = findNode(pt);
	if (n == NULL) {
		throw out_of_range("at");
	}
	return n->getVal();
}

template <size_t N, typename ElemType> const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
	Node<N, ElemType>* n = findNode(pt);
	if (n == NULL) {
		throw out_of_range("at");
	}
	return n->getVal();
}

template <size_t N, typename ElemType> void KDTree<N, ElemType>::search_subtree(priority_queue<pair<int, double>, vector<pair<int, double> >, CompareDist>& pq, const Node<N, ElemType>* curr, const Point<N>& key, int level, int k) const {
	if (curr == NULL) {
		return;
	}
	double dist = Distance(curr->getPoint(), key);
	// if pq has less than k elems in it, push curr on
	if (pq.size() < k) {
		pq.emplace((make_pair(curr->getVal(), dist)));
	}
	// otherwise, only push on if distance of curr is less than distance of max elem, and then pop max elem
	else if (dist < pq.top().second) {
		pq.pop();
		pq.emplace((make_pair(curr->getVal(), dist)));
	}
	bool left_child_searched = true;
	if (key[level] < curr->getPoint()[level]) {
		search_subtree(pq, curr->getLeftChild(), key, (level+1)%N, k);
	}
	else {
		search_subtree(pq, curr->getRightChild(), key, (level+1)%N, k);
		left_child_searched = false;
	}
	if ( (pq.size() < k) || (fabs(key[level] - curr->getPoint()[level]) < pq.top().second) ) {
		if (left_child_searched) {
			search_subtree(pq, curr->getRightChild(), key, (level+1)%N, k);
		}
		else {
			search_subtree(pq, curr->getLeftChild(), key, (level+1)%N, k);
		}
	}
}

template <size_t N, typename ElemType> void KDTree<N, ElemType>::kNN_Graph_recursive(int k, ElemType** neighbors, Node<N, ElemType>* node) const {
	if (node == NULL) { return; }
	priority_queue<pair<int, double>, vector<pair<int, double> >, CompareDist> pq;
	search_subtree(pq, root, node->getPoint(), 0, k);
	for (int i = k-1; i >= 0; --i) {
		neighbors[node->getVal()][i] = pq.top().first;
		pq.pop();
	}
	kNN_Graph_recursive(k, neighbors, node->getLeftChild());
	kNN_Graph_recursive(k, neighbors, node->getRightChild());
}

template <size_t N, typename ElemType> void KDTree<N, ElemType>::kNN_Graph(int k, ElemType** neighbors) const {
	kNN_Graph_recursive(k, neighbors, root);
	return;
}

#endif
