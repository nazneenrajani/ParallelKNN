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
#include "BoundedPQueue.h"
#include "Node.h"
#include <stdexcept>
#include <cmath>
using namespace std;

template <size_t N, typename ElemType> class KDTree {
public:
	/* Constructor: KDTree();
	 * Usage: KDTree myTree;
	 * ----------------------------------------------------
	 * Constructs an empty KDTree.
	 */
	KDTree(); // done

	/* Destructor: ~KDTree()
	 * Usage: (implicit)
	 * ----------------------------------------------------
	 * Cleans up all resources used by the KDTree.
	 */
	~KDTree(); // TODO

	/* KDTree(const KDTree& rhs);
	 * KDTree& operator= (const KDTree& rhs);
	 * Usage: KDTree one = two;
	 * Usage: one = two;
	 * -----------------------------------------------------
	 * Deep-copies the contents of another KDTree into this
	 * one.
	 */
//	KDTree(const KDTree& rhs); // needed?
//	KDTree& operator= (const KDTree& rhs); // needed?

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
	ElemType& at(const Point<N>& pt);
	const ElemType& at(const Point<N>& pt) const;

	/* ElemType kNNValue(const Point<N>& key, size_t k) const
	 * Usage: cout << kd.kNNValue(v, 3) << endl;
	 * ----------------------------------------------------
	 * Given a point v and an integer k, finds the k points
	 * in the KDTree nearest to v and returns the most common
	 * value associated with those points.  In the event of
	 * a tie, one of the most frequent value will be chosen.
	 */
	ElemType kNNValue(const Point<N>& key, size_t k) const;

private:
  Node<N, ElemType>* root;
  size_t n_items;
  Node<N, ElemType>* findNode(const Point<N>& pt) const; // returns a reference to the point if found; else returns NULL
};

/* * * * * Implementation Below This Point * * * * */

template <size_t N, typename ElemType> KDTree<N, ElemType>::KDTree() {
  n_items = 0;
  root = NULL;
}

template <size_t N, typename ElemType> KDTree<N, ElemType>::~KDTree() {
  // TODO: Fill this in.
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
//		return *d;
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


#endif
