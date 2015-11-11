//
//  Node.h
//  ttcr
//
//  Created by Bernard Giroux on 12-08-14.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//

#ifndef __NODE_H__
#define __NODE_H__

#include "ttRaisCourbes_t.h"

template<typename T>
class Node {
public:
	virtual T getTT(const int n) const = 0;
};

template<typename T>
class CompareNodePtr {
	
	// Overloaded operator for the priority queue, compare the "n"th traveltimes of two nodes.
private:
	int n;
public:
	CompareNodePtr(const int nn) : n(nn) {}
    bool operator()(const Node<T>* n1, const Node<T>* n2) const {
        //  The priority_queue must return the minimum time!!!
        return n1->getTT(n) > n2->getTT(n);
    }
};


#endif
