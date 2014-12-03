/*
 * Node.h
 *
 *  Created on: Nov 11, 2014
 *      Author: katemcardle
 */

#ifndef NODE_H_
#define NODE_H_

template <size_t N, typename ElemType> class Node {
public:
	Node(const Point<N>& pt, const ElemType& val);
	~Node();

	const Point<N>& getPoint() const;
	Node<N, ElemType>* getLeftChild() const;
	Node<N, ElemType>* getRightChild() const;
	const ElemType& getVal() const;
	ElemType& getVal();
	void setVal(const ElemType& val);
	bool add(const Point<N>& pt, const ElemType& value, int level);


private:
	ElemType val;
	Point<N> pt;
	Node* left;
	Node* right;
};

template <size_t N, typename ElemType> Node<N, ElemType>::Node(const Point<N>& pt, const ElemType& val) {
	this->val = val;
	this->pt = pt;
	left = NULL;
	right = NULL;
}

template <size_t N, typename ElemType> Node<N, ElemType>::~Node() {

}

template <size_t N, typename ElemType> const Point<N>& Node<N, ElemType>::getPoint() const {
	return pt;
}

template <size_t N, typename ElemType> ElemType& Node<N, ElemType>::getVal() {
//	printf("non-const getVal\n");
	return val;
}

template <size_t N, typename ElemType> const ElemType& Node<N, ElemType>::getVal() const {
//	printf("const getVal\n");
	return val;
}

template <size_t N, typename ElemType> void Node<N, ElemType>::setVal(const ElemType& val) {
	this->val = val;
}

template <size_t N, typename ElemType> Node<N, ElemType>* Node<N, ElemType>::getLeftChild() const {
	return left;
}

template <size_t N, typename ElemType> Node<N, ElemType>* Node<N, ElemType>::getRightChild() const {
	return right;
}

template <size_t N, typename ElemType> bool Node<N, ElemType>::add(const Point<N>& newPt, const ElemType& value, int level) {
	if (newPt == pt) {
		val = value;
		return false;
	}
	else if (newPt[level%N] < pt[level%N]) {
		if (left == NULL) {
			left = new Node<N, ElemType>(newPt, value);
			return true;
		}
		else {
			++level;
			return left->add(newPt, value, level);
		}
	}
	else {
		if (right == NULL) {
			right = new Node<N, ElemType>(newPt, value);
			return true;
		}
		else {
			++level;
			return right->add(newPt, value, level);
		}
	}
}

#endif /* NODE_H_ */
