//
//  Node2Dc.h
// ttRaisCourbes
//
// Created by Bernard Giroux on 08-04-24.
// Copyright 2008 Bernard Giroux.
//
//

//
// Copyright (C) 2008 Bernard Giroux.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __NODE2DC_H__
#define __NODE2DC_H__

#include <cmath>
#include <limits>

#include "Node.h"

template<typename T1, typename T2>
class Node2Dc : public Node<T1> {
public:
    Node2Dc(const int nt=1) :
	nThreads(nt),
    tt(0),
	x(0.0), z(0.0),
    gridIndex(std::numeric_limits<T2>::max()),
    nodeParent(0),
    cellParent(0),
    owners(0)
    {
		tt = new T1[nt];
		nodeParent = new T2[nt];
		cellParent = new T2[nt];

		for ( size_t n=0; n<nt; ++n ) {
			tt[n] = std::numeric_limits<T1>::max();
			nodeParent[n] = std::numeric_limits<T2>::max();
			cellParent[n] = std::numeric_limits<T2>::max();
		}
	}
	
    Node2Dc(const T1 xx, const T1 zz, const T2 index, const int nt=1) :
	nThreads(nt),
    tt(0),
	x(xx), z(zz),
    gridIndex(index),
    nodeParent(0),
    cellParent(0),
    owners(0)
    {
		tt = new T1[nt];
		nodeParent = new T2[nt];
		cellParent = new T2[nt];
		
		for ( size_t n=0; n<nt; ++n ) {
			tt[n] = std::numeric_limits<T1>::max();
			nodeParent[n] = std::numeric_limits<T2>::max();
			cellParent[n] = std::numeric_limits<T2>::max();
		}
	}
	
    
    Node2Dc(const T1 t, const T1 xx, const T1 zz, const int nt, const size_t i) :
	nThreads(nt),
    tt(0),
	x(xx), z(zz),
    gridIndex(std::numeric_limits<T2>::max()),
    nodeParent(0),
    cellParent(0),
    owners(std::vector<T2>(0))
    {
		tt = new T1[nt];
		nodeParent = new T2[nt];
		cellParent = new T2[nt];

		for ( size_t n=0; n<nt; ++n ) {
			tt[n] = std::numeric_limits<T1>::max();
			nodeParent[n] = std::numeric_limits<T2>::max();
			cellParent[n] = std::numeric_limits<T2>::max();
		}
		tt[i] = t;
	}
    
	Node2Dc(const Node2Dc<T1,T2>& node) :
	nThreads(node.nThreads), 
    tt(0),
	x(node.x), z(node.z),
    gridIndex(node.gridIndex),
    nodeParent(0),
    cellParent(0),
    owners(node.owners)
    {
		tt = new T1[nThreads];
		nodeParent = new T2[nThreads];
		cellParent = new T2[nThreads];
		
		for ( size_t n=0; n<nThreads; ++n ) {
			tt[n] = node.tt[n];
			nodeParent[n] = node.nodeParent[n];
			cellParent[n] = node.cellParent[n];
		}
	}
	

	~Node2Dc() {
		delete [] tt;
		delete [] nodeParent;
		delete [] cellParent;
	}
	
    void reinit(const size_t thread_no) { //=0) {
		tt[thread_no] = std::numeric_limits<T1>::max();
        nodeParent[thread_no] = std::numeric_limits<T2>::max();
        cellParent[thread_no] = std::numeric_limits<T2>::max();
    }
	
    T1 getTT(const int i) const { return tt[i]; }
    void setTT(const T1 t, const int i) { tt[i] = t; }

	void setXZindex(const T1 xx, const T1 zz, const T2 index) {
		x=xx; z=zz; gridIndex = index;  }

    T1 getX() const {
		return x;
	}
    void setX(const T1 xx) { x = xx; }
    
    T1 getZ() const { return z; }
    void setZ(const T1 zz) { z = zz; }
    
    T2 getGridIndex() const { return gridIndex; }
    void setGridIndex(const T2 index) { gridIndex = index; }
    
    T2 getNodeParent(const size_t i) const { return nodeParent[i]; }
    void setnodeParent(const T2 index, const size_t i) { nodeParent[i] = index; }
    
    T2 getCellParent(const size_t i) const { return cellParent[i]; }
    void setCellParent(const T2 index, const size_t i) { cellParent[i] = index; }
    
    void pushOwner(const T2 o) { owners.push_back(o); }    
    const std::vector<T2>& getOwners() const { return owners; }
    
    T1 getDistance( const Node2Dc<T1,T2>& node ) const {
        return sqrt( (x-node.x)*(x-node.x) + (z-node.z)*(z-node.z) );
    }
    
    T1 getDistance( const sxz<T1>& node ) const {
        return sqrt( (x-node.x)*(x-node.x) + (z-node.z)*(z-node.z) );
    }
    
    T1 getDistanceX( const sxz<T1>& node ) const {
        return fabs( x-node.x );
    }
    
    T1 getDistanceZ( const sxz<T1>& node ) const {
        return fabs( z-node.z );
    }
    
	// operator to test if same location
	bool operator==( const sxz<T1>& node ) const {
		return fabs(x-node.x)<small && fabs(z-node.z)<small;
	}
	
    size_t getSize() const {
        return sizeof(int) + nThreads*sizeof(T1) + 2*sizeof(T1) +
        (1+2*nThreads)*sizeof(T2) + owners.size() * sizeof(T2);
    }
	
	int getDimension() const { return 2; }

private:
	int nThreads;
	T1 *tt;                        // travel time
    T1 x;                          // x coordinate
    T1 z;                          // z coordinate
    T2 gridIndex;                  // index of this node in the list of the grid
    T2 *nodeParent;                // index of parent node of the ray
    T2 *cellParent;                // index of cell traversed by the ray
    std::vector<T2> owners;        // indices of cells touching the node
    
};


#endif
