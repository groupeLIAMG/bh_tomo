/*
 *  Grid2Dvti.h
 *  ttRaisCourbes
 *
 *  Created by Bernard Giroux on 13-03-06.
 *  Copyright 2013 Bernard Giroux.
 *
 *
 * 
 *
 */

/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __GRID2DVTI_H__
#define __GRID2DVTI_H__




#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>

#include "Node2Dc.h"




template<typename T1, typename T2>
class Grid2Dvti_sh {
 public:
	Grid2Dvti_sh(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
							 const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
							 const int nt=1);
	
	T1 getDx() const { return dx; }
	T1 getDz() const { return dz; }
	T1 getXmin() const { return xmin; }
	T1 getXmax() const { return xmax; }
	T1 getZmin() const { return zmin; }
	T1 getZmax() const { return zmax; }
	T2 getNcellx() const { return nCellx; }
	T2 getNcellz() const { return nCellz; }
	T2 getNsnx() const { return nsnx; }
	T2 getNsnz() const { return nsnz; }
    
	T2 getNumberOfCells() const { return Vs0.size(); }
    
	virtual void setVs0(const T1 v) {
		for ( size_t n=0; n<Vs0.size(); ++n ) {
			Vs0[n] = v;
		}
	}
    
	virtual int setVs0(const T1 *v, const size_t nv) {
		if ( Vs0.size() != nv ) {
			std::cerr << "Error:  vectors of incompatible size.";
			return 1;
		}
		for ( size_t n=0; n<Vs0.size(); ++n ) {
			Vs0[n] = v[n];
		}
		return 0;
	}
    
	virtual void setGamma(const T1 g) {
		for ( size_t n=0; n<gamma.size(); ++n ) {
			gamma[n] = g;
		}
	}
    
	virtual int setGamma(const T1 *g, const size_t ng) {
		if ( gamma.size() != ng ) {
			std::cerr << "Error:  vectors of incompatible size.";
			return 1;
		}
		for ( size_t n=0; n<gamma.size(); ++n ) {
			gamma[n] = g[n];
		}
		return 0;
	}
    
	int raytrace(const std::vector<sxz<T1> >& Tx,
							 const std::vector<T1>& t0, 
							 const std::vector<sxz<T1> >& Rx,
							 std::vector<T1>& traveltimes,
							 const int threadNo=0) const;
    
	int raytrace(const std::vector<sxz<T1> >& Tx,
							 const std::vector<T1>& t0, 
							 const std::vector<sxz<T1> >& Rx,
							 std::vector<T1>& traveltimes,
							 std::vector<std::vector<sxz<double> > >& r_data,
							 const int threadNo=0) const;
    
    
	T2 getCellNo(const sxz<T1>& pt) const {
		T2 nx = static_cast<T2>( small + (pt.x-xmin)/dx );
		T2 nz = static_cast<T2>( small + (pt.z-zmin)/dz );
		return nx*nCellz + nz;
	}

   
    
 protected:
	int nThreads;
	T1 dx;           // cell size in x
	T1 dz;           // cell size in z
	T1 xmin;         // x origin of the grid
	T1 zmin;         // z origin of the grid
	T1 xmax;         // x end of the grid
	T1 zmax;         // z end of the grid
	const T1 small;
	T2 nCellx;  // number of cells in x
	T2 nCellz;  // number of cells in x
	T2 nsnx;    // number of secondary nodes in x
	T2 nsnz;    // number of secondary nodes in z
	T2 nsgx;    // number of subgrid cells in x
	T2 nsgz;    // number of subgrid cells in z
    
	mutable std::vector<Node2Dc<T1,T2> > nodes;
    
	std::vector<T1> Vs0;   // column-wise (z axis) velocity vector of the cells
	std::vector<T1> gamma; // column-wise (z axis) anisotropy parameter vector of the cells
	std::vector<std::vector<T2> > neighbors;  // nodes common to a cell
    
	void buildGridNodes();
	void buildGridNeighbors();
	
	void propagate(std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
								 CompareNodePtr<T1> >& queue,
								 std::vector<bool>& inQueue,
								 std::vector<bool>& frozen,
								 const int threadNo) const;
	
	virtual T1 computeDt(const Node2Dc<T1,T2>& source, const sxz<T1>& node,
											 const T2 cellNo) const {
		// theta: angle w/r to vertical axis
		T1 theta = atan2(node.x - source.getX(), node.z - source.getZ());
		T1 v = Vs0[cellNo] * sqrt(1. + 2.*gamma[cellNo]*sin(theta)*sin(theta));
		return source.getDistance( node ) / v;
	}
    
	virtual T1 computeDt(const Node2Dc<T1,T2>& source, const Node2Dc<T1,T2>& node,
											 const T2 cellNo) const {
		// theta: angle w/r to vertical axis
		T1 theta = atan2(node.getX() - source.getX(), node.getZ() - source.getZ());
		T1 v = Vs0[cellNo] * sqrt(1. + 2.*gamma[cellNo]*sin(theta)*sin(theta));
		return source.getDistance( node ) / v;
	}
    
	int check_pts(const std::vector<sxz<T1> >&) const;
    
	void initQueue(const std::vector<sxz<T1> >& Tx,
								 const std::vector<T1>& t0, 
								 std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
								 CompareNodePtr<T1> >& queue,
								 std::vector<Node2Dc<T1,T2> >& txNodes,
								 std::vector<bool>& inQueue,
								 std::vector<bool>& frozen,
								 const int threadNo) const;
    
   
	T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Dc<T1,T2> >& nodes,
									 const int threadNo) const;
    
	T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Dc<T1,T2> >& nodes,
									 T2& nodeParentRx, T2& cellParentRx,
									 const int threadNo) const;

	void save(const char filename[]) const;
    
	
 private:
	Grid2Dvti_sh() {}
	Grid2Dvti_sh(const Grid2Dvti_sh<T1,T2>& g) {}
	Grid2Dvti_sh<T1,T2>& operator=(const Grid2Dvti_sh<T1,T2>& g) {}
    
};

template<typename T1, typename T2>
	Grid2Dvti_sh<T1,T2>::Grid2Dvti_sh(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
																		const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
																		const int nt) : nThreads(nt),
	dx(ddx), dz(ddz), xmin(minx), zmin(minz), xmax(minx+nx*ddx), zmax(minz+nz*ddz),
	small(1.0e-10), nCellx(nx), nCellz(nz), nsnx(nnx), nsnz(nnz), nsgx(0), nsgz(0),
	nodes(std::vector<Node2Dc<T1,T2> >( // noeuds secondaires
																		 nCellx*nsnx*(nCellz+1) + nCellz*nsnz*(nCellx+1) +
																		 // noeuds primaires
																		 (nCellx+1) * (nCellz+1), Node2Dc<T1,T2>(nt) )),
	Vs0(std::vector<T1>(nCellx*nCellz)), gamma(std::vector<T1>(nCellx*nCellz)), 
	neighbors(std::vector<std::vector<T2> >(nCellx*nCellz))
																		{
																			buildGridNodes();
																			buildGridNeighbors();
																		}

template<typename T1, typename T2>
	void Grid2Dvti_sh<T1,T2>::buildGridNodes() {
	T1 dxs = dx/(nsnx+1);
	T1 dzs = dz/(nsnz+1);
    
	T2 cell_upLeft = std::numeric_limits<T2>::max();
	T2 cell_upRight = std::numeric_limits<T2>::max();
	T2 cell_downLeft = 0;
	T2 cell_downRight = 0;
    
	for ( T2 n=0, nc=0; nc<=nCellx; ++nc ) {
        
		double x = xmin + nc*dx;
        
		for ( T2 nr=0; nr<=nCellz; ++nr ) {
            
			double z = zmin + nr*dz;
            
			if ( nr < nCellz && nc < nCellx ) {
				cell_downRight = nc*nCellz + nr;
			}
			else {
				cell_downRight = std::numeric_limits<T2>::max();
			}
            
			if ( nr > 0 && nc < nCellx ) {
				cell_upRight = nc*nCellz + nr - 1;
			}
			else {
				cell_upRight = std::numeric_limits<T2>::max();
			}
            
			if ( nr < nCellz && nc > 0 ) {
				cell_downLeft = (nc-1)*nCellz + nr;
			}
			else {
				cell_downLeft = std::numeric_limits<T2>::max();
			}
            
			if ( nr > 0 && nc > 0 ) {
				cell_upLeft = (nc-1)*nCellz + nr - 1;
			}
			else {
				cell_upLeft = std::numeric_limits<T2>::max();
			}
            
			if ( cell_upLeft != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_upLeft );
			}
			if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_downLeft );
			}
			if ( cell_upRight != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_upRight );
			}
			if ( cell_downRight != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_downRight );
			}
            
			nodes[n].setX( x );
			nodes[n].setZ( z );
			nodes[n].setGridIndex( n );
            
			++n;
            
			// secondary nodes on the vertical
			if ( nr < nCellz ) {
				for (T2 ns=0; ns<nsnz; ++ns, ++n ) {
                    
					double zsv = zmin + nr*dz + (ns+1)*dzs;
                    
					if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_downLeft );
					}
					if ( cell_downRight != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_downRight );
					}
                    
					nodes[n].setX( x );
					nodes[n].setZ( zsv );
					nodes[n].setGridIndex( n );
				}
			}
            
			// secondary nodes on the horizontal
			if ( nc < nCellx ) {
				for ( T2 ns=0; ns<nsnx; ++ns, ++n ) {
                    
					double xsh = xmin + nc*dx + (ns+1)*dxs;
                    
					if ( cell_upRight != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_upRight );
					}
					if ( cell_downRight != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_downRight );
					}
                    
					nodes[n].setX( xsh );
					nodes[n].setZ( z );
					nodes[n].setGridIndex( n );
				}
			}
		}
	}
}

template<typename T1, typename T2>
	void Grid2Dvti_sh<T1,T2>::buildGridNeighbors() {
	for ( T2 n=0; n<nodes.size(); ++n ) {
		for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
			neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
		}
	}
}


template<typename T1, typename T2>
	int Grid2Dvti_sh<T1,T2>::check_pts(const std::vector<sxz<T1> >& pts) const {
	for (size_t n=0; n<pts.size(); ++n) {
		if ( pts[n].x < xmin || pts[n].x > xmax ||
				 pts[n].z < zmin || pts[n].z > zmax ) {
			std::cerr << "Error: point no " << (n+1)
								<< " outside the grid.\n";
			return 1;
		}
	}
	return 0;
}


template<typename T1, typename T2>
	void Grid2Dvti_sh<T1,T2>::initQueue(const std::vector<sxz<T1> >& Tx,
																			const std::vector<T1>& t0, 
																			std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
																			CompareNodePtr<T1> >& queue,
																			std::vector<Node2Dc<T1,T2> >& txNodes,
																			std::vector<bool>& inQueue,
																			std::vector<bool>& frozen,
																			const int threadNo) const {
    
	for (size_t n=0; n<Tx.size(); ++n) {
		bool found = false;
		for ( size_t nn=0; nn<nodes.size(); ++nn ) {
			if ( nodes[nn] == Tx[n] ) {
				found = true;
				nodes[nn].setTT( t0[n], threadNo );
				queue.push( &(nodes[nn]) );
				inQueue[nn] = true;
				frozen[nn] = true;
				break;
			}
		}
		if ( found==false ) {
			txNodes.push_back( Node2Dc<T1,T2>(t0[n], Tx[n].x, Tx[n].z, nThreads, threadNo) );
			// we belong to cell index no
			txNodes.back().pushOwner( getCellNo(Tx[n]) );
			txNodes.back().setGridIndex( nodes.size()+txNodes.size()-1 );
            
			queue.push( &(txNodes.back()) );
			inQueue.push_back( true );
			frozen.push_back( true );
		}
	}
}


template<typename T1, typename T2>
	int Grid2Dvti_sh<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
																		const std::vector<T1>& t0, 
																		const std::vector<sxz<T1> >& Rx,
																		std::vector<T1>& traveltimes,
																		const int threadNo) const {
    
	//	std::cout << "   running in thread no " << threadNo << std::endl;
	if ( check_pts(Tx) == 1 ) return 1;
	if ( check_pts(Rx) == 1 ) return 1;
    
	for ( size_t n=0; n<nodes.size(); ++n ) {
		nodes[n].reinit( threadNo );
	}
    
	CompareNodePtr<T1> cmp(threadNo);
	std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    
	std::vector<Node2Dc<T1,T2> > txNodes;
	std::vector<bool> inQueue( nodes.size(), false );
	std::vector<bool> frozen( nodes.size(), false );
    
	initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
	propagate(queue, inQueue, frozen, threadNo);
    
	if ( traveltimes.size() != Rx.size() ) {
		traveltimes.resize( Rx.size() );
	}
    
	for (size_t n=0; n<Rx.size(); ++n) {
		traveltimes[n] = getTraveltime(Rx[n], nodes, threadNo);
	}
	return 0;
}

template<typename T1, typename T2>
	int Grid2Dvti_sh<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
																		const std::vector<T1>& t0, 
																		const std::vector<sxz<T1> >& Rx,
																		std::vector<T1>& traveltimes,
																		std::vector<std::vector<sxz<double> > >& r_data,
																		const int threadNo) const {
    
	if ( check_pts(Tx) == 1 ) return 1;
	if ( check_pts(Rx) == 1 ) return 1;
    
	for ( size_t n=0; n<nodes.size(); ++n ) {
		nodes[n].reinit( threadNo );
	}
    
	CompareNodePtr<T1> cmp(threadNo);
	std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
	std::vector<Node2Dc<T1,T2> > txNodes;
	std::vector<bool> inQueue( nodes.size(), false );
	std::vector<bool> frozen( nodes.size(), false );
    
	initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
	propagate(queue, inQueue, frozen, threadNo);
    
	if ( traveltimes.size() != Rx.size() ) {
		traveltimes.resize( Rx.size() );
	}
	if ( r_data.size() != Rx.size() ) {
		r_data.resize( Rx.size() );
	}
	for ( size_t ni=0; ni<r_data.size(); ++ni ) {
		r_data[ni].resize( 0 );
	}
	T2 nodeParentRx;
	T2 cellParentRx;
    
	for (size_t n=0; n<Rx.size(); ++n) {
        
		traveltimes[n] = getTraveltime(Rx[n], nodes, nodeParentRx, cellParentRx,
																	 threadNo);
        
		// Rx are in nodes (not txNodes)
		std::vector<Node2Dc<T1,T2> > *node_p;
		node_p = &nodes;
        
		std::vector<sxz<double> > r_tmp;
		T2 iChild, iParent = nodeParentRx;
		sxz<double> child;
		
		// store the son's coord 
		child.x = Rx[n].x;
		child.z = Rx[n].z;

		while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
 			
			r_tmp.push_back( child );

			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
			
			// grand'pa is now papa
			iParent = (*node_p)[iChild].getNodeParent(threadNo);
			if ( iParent >= nodes.size() ) {
				node_p = &txNodes;
				iParent -= nodes.size();
			}
			else {
				node_p = &nodes;
			}
		}
		
		// parent is now at Tx
		r_tmp.push_back( child );
			
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );
		
       
		// the order should be from Tx to Rx, so we reorder...
		iParent = r_tmp.size();
		r_data[n].resize( r_tmp.size() );
		for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
			r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
			r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
		}
	}
	return 0;
}





template<typename T1, typename T2>
	void Grid2Dvti_sh<T1,T2>::propagate( std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
																			 CompareNodePtr<T1> >& queue,
																			 std::vector<bool>& inQueue,
																			 std::vector<bool>& frozen,
																			 const int threadNo) const {
    
	while ( !queue.empty() ) {
		const Node2Dc<T1,T2>* source = queue.top();
		queue.pop();
		inQueue[ source->getGridIndex() ] = false;
        
		for ( size_t no=0; no<source->getOwners().size(); ++no ) {
            
			T2 cellNo = source->getOwners()[no];
            
			for ( size_t k=0; k< neighbors[cellNo].size(); ++k ) {
				T2 neibNo = neighbors[cellNo][k];
				if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
					continue;
				}
                
				// compute dt
				T1 dt = computeDt(*source, nodes[neibNo], cellNo);
				
				if ( source->getTT(threadNo)+dt < nodes[neibNo].getTT(threadNo) ) {
					nodes[neibNo].setTT( source->getTT(threadNo)+dt, threadNo );
					nodes[neibNo].setnodeParent( source->getGridIndex(), threadNo );
					nodes[neibNo].setCellParent( cellNo, threadNo );
                    
					if ( !inQueue[neibNo] ) {
						queue.push( &(nodes[neibNo]) );
						inQueue[neibNo] = true;
					}
				}
			}
		}
	}
}





template<typename T1, typename T2>
	T1 Grid2Dvti_sh<T1,T2>::getTraveltime(const sxz<T1>& Rx,
																				const std::vector<Node2Dc<T1,T2> >& nodes,
																				const int threadNo) const {
    
	for ( size_t nn=0; nn<nodes.size(); ++nn ) {
		if ( nodes[nn] == Rx ) {
			return nodes[nn].getTT(threadNo);
		}
	}
    
	T2 cellNo = getCellNo( Rx );
	T2 neibNo = neighbors[cellNo][0];
	T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
	T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
	for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
		neibNo = neighbors[cellNo][k];
		dt = computeDt(nodes[neibNo], Rx, cellNo);
		if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
			traveltime =  nodes[neibNo].getTT(threadNo)+dt;
		}
	}
	return traveltime;
}


template<typename T1, typename T2>
	T1 Grid2Dvti_sh<T1,T2>::getTraveltime(const sxz<T1>& Rx,
																				const std::vector<Node2Dc<T1,T2> >& nodes,
																				T2& nodeParentRx, T2& cellParentRx,
																				const int threadNo) const {
    
	for ( size_t nn=0; nn<nodes.size(); ++nn ) {
		if ( nodes[nn] == Rx ) {
			nodeParentRx = nodes[nn].getNodeParent(threadNo);
			cellParentRx = nodes[nn].getCellParent(threadNo);
			return nodes[nn].getTT(threadNo);
		}
	}
    
	T2 cellNo = getCellNo( Rx );
	T2 neibNo = neighbors[cellNo][0];
	T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
	T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
	nodeParentRx = neibNo;
	cellParentRx = cellNo;
	for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
		neibNo = neighbors[cellNo][k];
		dt = computeDt(nodes[neibNo], Rx, cellNo);
		if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
			traveltime =  nodes[neibNo].getTT(threadNo)+dt;
			nodeParentRx = neibNo;
		}
	}
	return traveltime;
}









template<typename T1, typename T2>
class Grid2Dvti_psv {
 public:
	Grid2Dvti_psv(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
							 const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
							 const int nt=1);
	
	T1 getDx() const { return dx; }
	T1 getDz() const { return dz; }
	T1 getXmin() const { return xmin; }
	T1 getXmax() const { return xmax; }
	T1 getZmin() const { return zmin; }
	T1 getZmax() const { return zmax; }
	T2 getNcellx() const { return nCellx; }
	T2 getNcellz() const { return nCellz; }
	T2 getNsnx() const { return nsnx; }
	T2 getNsnz() const { return nsnz; }
    
	T2 getNumberOfCells() const { return Vs0.size(); }
    
	virtual int setVp0(const T1 *v, const size_t nv) {
		if ( Vp0.size() != nv ) {
			std::cerr << "Error:  vectors of incompatible size.";
			return 1;
		}
		for ( size_t n=0; n<Vp0.size(); ++n ) {
			Vp0[n] = v[n];
		}
		return 0;
	}
    
	virtual int setVs0(const T1 *v, const size_t nv) {
		if ( Vs0.size() != nv ) {
			std::cerr << "Error:  vectors of incompatible size.";
			return 1;
		}
		for ( size_t n=0; n<Vs0.size(); ++n ) {
			Vs0[n] = v[n];
		}
		return 0;
	}
    
	virtual int setEpsilon(const T1 *g, const size_t ng) {
		if ( epsilon.size() != ng ) {
			std::cerr << "Error:  vectors of incompatible size.";
			return 1;
		}
		for ( size_t n=0; n<epsilon.size(); ++n ) {
			epsilon[n] = g[n];
		}
		return 0;
	}
	
	virtual int setDelta(const T1 *g, const size_t ng) {
		if ( delta.size() != ng ) {
			std::cerr << "Error:  vectors of incompatible size.";
			return 1;
		}
		for ( size_t n=0; n<delta.size(); ++n ) {
			delta[n] = g[n];
		}
		return 0;
	}

	void setPhase(const int p) {
		if ( p==1 ) sign = 1.;  // P wave
		else sign = -1.;        // SV wave
	}

	int raytrace(const std::vector<sxz<T1> >& Tx,
							 const std::vector<T1>& t0, 
							 const std::vector<sxz<T1> >& Rx,
							 std::vector<T1>& traveltimes,
							 const int threadNo=0) const;
    
	int raytrace(const std::vector<sxz<T1> >& Tx,
							 const std::vector<T1>& t0, 
							 const std::vector<sxz<T1> >& Rx,
							 std::vector<T1>& traveltimes,
							 std::vector<std::vector<sxz<double> > >& r_data,
							 const int threadNo=0) const;
    
    
	T2 getCellNo(const sxz<T1>& pt) const {
		T2 nx = static_cast<T2>( small + (pt.x-xmin)/dx );
		T2 nz = static_cast<T2>( small + (pt.z-zmin)/dz );
		return nx*nCellz + nz;
	}

   
    
 protected:
	int nThreads;
	T1 sign;
	T1 dx;           // cell size in x
	T1 dz;           // cell size in z
	T1 xmin;         // x origin of the grid
	T1 zmin;         // z origin of the grid
	T1 xmax;         // x end of the grid
	T1 zmax;         // z end of the grid
	const T1 small;
	T2 nCellx;  // number of cells in x
	T2 nCellz;  // number of cells in x
	T2 nsnx;    // number of secondary nodes in x
	T2 nsnz;    // number of secondary nodes in z
	T2 nsgx;    // number of subgrid cells in x
	T2 nsgz;    // number of subgrid cells in z
    
	mutable std::vector<Node2Dc<T1,T2> > nodes;
    
	std::vector<T1> Vp0;   // column-wise (z axis) velocity vector of the cells
	std::vector<T1> Vs0;   // column-wise (z axis) velocity vector of the cells
	std::vector<T1> epsilon; // column-wise (z axis) anisotropy parameter vector of the cells
	std::vector<T1> delta;   // column-wise (z axis) anisotropy parameter vector of the cells
	std::vector<std::vector<T2> > neighbors;  // nodes common to a cell
    
	void buildGridNodes();
	void buildGridNeighbors();
	
	void propagate(std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
								 CompareNodePtr<T1> >& queue,
								 std::vector<bool>& inQueue,
								 std::vector<bool>& frozen,
								 const int threadNo) const;
	
	virtual T1 computeDt(const Node2Dc<T1,T2>& source, const sxz<T1>& node,
											 const T2 cellNo) const {
		// theta: angle w/r to vertical axis
		T1 theta = atan2(node.x - source.getX(), node.z - source.getZ());
		T1 f = 1. - (Vs0[cellNo]*Vs0[cellNo]) / (Vp0[cellNo]*Vp0[cellNo]);

		T1 tmp = 1. + (2.*epsilon[cellNo]*sin(theta)*sin(theta)) / f;

		tmp = 1. + epsilon[cellNo]*sin(theta)*sin(theta) - f/2. +
			sign*f/2.*sqrt( tmp*tmp - (2.*(epsilon[cellNo]-delta[cellNo])*sin(2.*theta)*sin(2.*theta))/f );

		T1 v = Vp0[cellNo] * sqrt( tmp );
		return source.getDistance( node ) / v;
	}
    
	virtual T1 computeDt(const Node2Dc<T1,T2>& source, const Node2Dc<T1,T2>& node,
											 const T2 cellNo) const {
		// theta: angle w/r to vertical axis
		T1 theta = atan2(node.getX() - source.getX(), node.getZ() - source.getZ());
		T1 f = 1. - (Vs0[cellNo]*Vs0[cellNo]) / (Vp0[cellNo]*Vp0[cellNo]);

		T1 tmp = 1. + (2.*epsilon[cellNo]*sin(theta)*sin(theta)) / f;

		tmp = 1. + epsilon[cellNo]*sin(theta)*sin(theta) - f/2. +
			sign*f/2.*sqrt( tmp*tmp - (2.*(epsilon[cellNo]-delta[cellNo])*sin(2.*theta)*sin(2.*theta))/f );

		T1 v = Vp0[cellNo] * sqrt( tmp );
		return source.getDistance( node ) / v;
	}
	
	int check_pts(const std::vector<sxz<T1> >&) const;
    
	void initQueue(const std::vector<sxz<T1> >& Tx,
								 const std::vector<T1>& t0, 
								 std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
								 CompareNodePtr<T1> >& queue,
								 std::vector<Node2Dc<T1,T2> >& txNodes,
								 std::vector<bool>& inQueue,
								 std::vector<bool>& frozen,
								 const int threadNo) const;
    
   
	T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Dc<T1,T2> >& nodes,
									 const int threadNo) const;
    
	T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Dc<T1,T2> >& nodes,
									 T2& nodeParentRx, T2& cellParentRx,
									 const int threadNo) const;

	void save(const char filename[]) const;
    
	
 private:
	Grid2Dvti_psv() {}
	Grid2Dvti_psv(const Grid2Dvti_psv<T1,T2>& g) {}
	Grid2Dvti_psv<T1,T2>& operator=(const Grid2Dvti_psv<T1,T2>& g) {}
    
};

template<typename T1, typename T2>
	Grid2Dvti_psv<T1,T2>::Grid2Dvti_psv(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
																		const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
																			const int nt) : nThreads(nt), sign(1.), 
	dx(ddx), dz(ddz), xmin(minx), zmin(minz), xmax(minx+nx*ddx), zmax(minz+nz*ddz),
	small(1.0e-10), nCellx(nx), nCellz(nz), nsnx(nnx), nsnz(nnz), nsgx(0), nsgz(0),
	nodes(std::vector<Node2Dc<T1,T2> >( // noeuds secondaires
																		 nCellx*nsnx*(nCellz+1) + nCellz*nsnz*(nCellx+1) +
																		 // noeuds primaires
																		 (nCellx+1) * (nCellz+1), Node2Dc<T1,T2>(nt) )),
	Vp0(std::vector<T1>(nCellx*nCellz)), Vs0(std::vector<T1>(nCellx*nCellz)), 
	epsilon(std::vector<T1>(nCellx*nCellz)), delta(std::vector<T1>(nCellx*nCellz)),
	neighbors(std::vector<std::vector<T2> >(nCellx*nCellz))
																		{
																			buildGridNodes();
																			buildGridNeighbors();
																		}

template<typename T1, typename T2>
	void Grid2Dvti_psv<T1,T2>::buildGridNodes() {
	T1 dxs = dx/(nsnx+1);
	T1 dzs = dz/(nsnz+1);
    
	T2 cell_upLeft = std::numeric_limits<T2>::max();
	T2 cell_upRight = std::numeric_limits<T2>::max();
	T2 cell_downLeft = 0;
	T2 cell_downRight = 0;
    
	for ( T2 n=0, nc=0; nc<=nCellx; ++nc ) {
        
		double x = xmin + nc*dx;
        
		for ( T2 nr=0; nr<=nCellz; ++nr ) {
            
			double z = zmin + nr*dz;
            
			if ( nr < nCellz && nc < nCellx ) {
				cell_downRight = nc*nCellz + nr;
			}
			else {
				cell_downRight = std::numeric_limits<T2>::max();
			}
            
			if ( nr > 0 && nc < nCellx ) {
				cell_upRight = nc*nCellz + nr - 1;
			}
			else {
				cell_upRight = std::numeric_limits<T2>::max();
			}
            
			if ( nr < nCellz && nc > 0 ) {
				cell_downLeft = (nc-1)*nCellz + nr;
			}
			else {
				cell_downLeft = std::numeric_limits<T2>::max();
			}
            
			if ( nr > 0 && nc > 0 ) {
				cell_upLeft = (nc-1)*nCellz + nr - 1;
			}
			else {
				cell_upLeft = std::numeric_limits<T2>::max();
			}
            
			if ( cell_upLeft != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_upLeft );
			}
			if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_downLeft );
			}
			if ( cell_upRight != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_upRight );
			}
			if ( cell_downRight != std::numeric_limits<T2>::max() ) {
				nodes[n].pushOwner( cell_downRight );
			}
            
			nodes[n].setX( x );
			nodes[n].setZ( z );
			nodes[n].setGridIndex( n );
            
			++n;
            
			// secondary nodes on the vertical
			if ( nr < nCellz ) {
				for (T2 ns=0; ns<nsnz; ++ns, ++n ) {
                    
					double zsv = zmin + nr*dz + (ns+1)*dzs;
                    
					if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_downLeft );
					}
					if ( cell_downRight != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_downRight );
					}
                    
					nodes[n].setX( x );
					nodes[n].setZ( zsv );
					nodes[n].setGridIndex( n );
				}
			}
            
			// secondary nodes on the horizontal
			if ( nc < nCellx ) {
				for ( T2 ns=0; ns<nsnx; ++ns, ++n ) {
                    
					double xsh = xmin + nc*dx + (ns+1)*dxs;
                    
					if ( cell_upRight != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_upRight );
					}
					if ( cell_downRight != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_downRight );
					}
                    
					nodes[n].setX( xsh );
					nodes[n].setZ( z );
					nodes[n].setGridIndex( n );
				}
			}
		}
	}
}

template<typename T1, typename T2>
	void Grid2Dvti_psv<T1,T2>::buildGridNeighbors() {
	for ( T2 n=0; n<nodes.size(); ++n ) {
		for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
			neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
		}
	}
}


template<typename T1, typename T2>
	int Grid2Dvti_psv<T1,T2>::check_pts(const std::vector<sxz<T1> >& pts) const {
	for (size_t n=0; n<pts.size(); ++n) {
		if ( pts[n].x < xmin || pts[n].x > xmax ||
				 pts[n].z < zmin || pts[n].z > zmax ) {
			std::cerr << "Error: point no " << (n+1)
								<< " outside the grid.\n";
			return 1;
		}
	}
	return 0;
}


template<typename T1, typename T2>
	void Grid2Dvti_psv<T1,T2>::initQueue(const std::vector<sxz<T1> >& Tx,
																			const std::vector<T1>& t0, 
																			std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
																			CompareNodePtr<T1> >& queue,
																			std::vector<Node2Dc<T1,T2> >& txNodes,
																			std::vector<bool>& inQueue,
																			std::vector<bool>& frozen,
																			const int threadNo) const {
    
	for (size_t n=0; n<Tx.size(); ++n) {
		bool found = false;
		for ( size_t nn=0; nn<nodes.size(); ++nn ) {
			if ( nodes[nn] == Tx[n] ) {
				found = true;
				nodes[nn].setTT( t0[n], threadNo );
				queue.push( &(nodes[nn]) );
				inQueue[nn] = true;
				frozen[nn] = true;
				break;
			}
		}
		if ( found==false ) {
			txNodes.push_back( Node2Dc<T1,T2>(t0[n], Tx[n].x, Tx[n].z, nThreads, threadNo) );
			// we belong to cell index no
			txNodes.back().pushOwner( getCellNo(Tx[n]) );
			txNodes.back().setGridIndex( nodes.size()+txNodes.size()-1 );
            
			queue.push( &(txNodes.back()) );
			inQueue.push_back( true );
			frozen.push_back( true );
		}
	}
}


template<typename T1, typename T2>
	int Grid2Dvti_psv<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
																		const std::vector<T1>& t0, 
																		const std::vector<sxz<T1> >& Rx,
																		std::vector<T1>& traveltimes,
																		const int threadNo) const {
    
	//	std::cout << "   running in thread no " << threadNo << std::endl;
	if ( check_pts(Tx) == 1 ) return 1;
	if ( check_pts(Rx) == 1 ) return 1;
    
	for ( size_t n=0; n<nodes.size(); ++n ) {
		nodes[n].reinit( threadNo );
	}
    
	CompareNodePtr<T1> cmp(threadNo);
	std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    
	std::vector<Node2Dc<T1,T2> > txNodes;
	std::vector<bool> inQueue( nodes.size(), false );
	std::vector<bool> frozen( nodes.size(), false );
    
	initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
	propagate(queue, inQueue, frozen, threadNo);
    
	if ( traveltimes.size() != Rx.size() ) {
		traveltimes.resize( Rx.size() );
	}
    
	for (size_t n=0; n<Rx.size(); ++n) {
		traveltimes[n] = getTraveltime(Rx[n], nodes, threadNo);
	}
	return 0;
}

template<typename T1, typename T2>
	int Grid2Dvti_psv<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
																		const std::vector<T1>& t0, 
																		const std::vector<sxz<T1> >& Rx,
																		std::vector<T1>& traveltimes,
																		std::vector<std::vector<sxz<double> > >& r_data,
																		const int threadNo) const {
    
	if ( check_pts(Tx) == 1 ) return 1;
	if ( check_pts(Rx) == 1 ) return 1;
    
	for ( size_t n=0; n<nodes.size(); ++n ) {
		nodes[n].reinit( threadNo );
	}
    
	CompareNodePtr<T1> cmp(threadNo);
	std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
	std::vector<Node2Dc<T1,T2> > txNodes;
	std::vector<bool> inQueue( nodes.size(), false );
	std::vector<bool> frozen( nodes.size(), false );
    
	initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
	propagate(queue, inQueue, frozen, threadNo);
    
	if ( traveltimes.size() != Rx.size() ) {
		traveltimes.resize( Rx.size() );
	}
	if ( r_data.size() != Rx.size() ) {
		r_data.resize( Rx.size() );
	}
	for ( size_t ni=0; ni<r_data.size(); ++ni ) {
		r_data[ni].resize( 0 );
	}
	T2 nodeParentRx;
	T2 cellParentRx;
    
	for (size_t n=0; n<Rx.size(); ++n) {
        
		traveltimes[n] = getTraveltime(Rx[n], nodes, nodeParentRx, cellParentRx,
																	 threadNo);
        
		// Rx are in nodes (not txNodes)
		std::vector<Node2Dc<T1,T2> > *node_p;
		node_p = &nodes;
        
		std::vector<sxz<double> > r_tmp;
		T2 iChild, iParent = nodeParentRx;
		sxz<double> child;
		
		// store the son's coord 
		child.x = Rx[n].x;
		child.z = Rx[n].z;

		while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
 			
			r_tmp.push_back( child );

			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
			
			// grand'pa is now papa
			iParent = (*node_p)[iChild].getNodeParent(threadNo);
			if ( iParent >= nodes.size() ) {
				node_p = &txNodes;
				iParent -= nodes.size();
			}
			else {
				node_p = &nodes;
			}
		}
		
		// parent is now at Tx
		r_tmp.push_back( child );
			
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );
		
       
		// the order should be from Tx to Rx, so we reorder...
		iParent = r_tmp.size();
		r_data[n].resize( r_tmp.size() );
		for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
			r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
			r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
		}
	}
	return 0;
}





template<typename T1, typename T2>
	void Grid2Dvti_psv<T1,T2>::propagate( std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
																			 CompareNodePtr<T1> >& queue,
																			 std::vector<bool>& inQueue,
																			 std::vector<bool>& frozen,
																			 const int threadNo) const {
    
	while ( !queue.empty() ) {
		const Node2Dc<T1,T2>* source = queue.top();
		queue.pop();
		inQueue[ source->getGridIndex() ] = false;
        
		for ( size_t no=0; no<source->getOwners().size(); ++no ) {
            
			T2 cellNo = source->getOwners()[no];
            
			for ( size_t k=0; k< neighbors[cellNo].size(); ++k ) {
				T2 neibNo = neighbors[cellNo][k];
				if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
					continue;
				}
                
				// compute dt
				T1 dt = computeDt(*source, nodes[neibNo], cellNo);
				
				if ( source->getTT(threadNo)+dt < nodes[neibNo].getTT(threadNo) ) {
					nodes[neibNo].setTT( source->getTT(threadNo)+dt, threadNo );
					nodes[neibNo].setnodeParent( source->getGridIndex(), threadNo );
					nodes[neibNo].setCellParent( cellNo, threadNo );
                    
					if ( !inQueue[neibNo] ) {
						queue.push( &(nodes[neibNo]) );
						inQueue[neibNo] = true;
					}
				}
			}
		}
	}
}





template<typename T1, typename T2>
	T1 Grid2Dvti_psv<T1,T2>::getTraveltime(const sxz<T1>& Rx,
																				const std::vector<Node2Dc<T1,T2> >& nodes,
																				const int threadNo) const {
    
	for ( size_t nn=0; nn<nodes.size(); ++nn ) {
		if ( nodes[nn] == Rx ) {
			return nodes[nn].getTT(threadNo);
		}
	}
    
	T2 cellNo = getCellNo( Rx );
	T2 neibNo = neighbors[cellNo][0];
	T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
	T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
	for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
		neibNo = neighbors[cellNo][k];
		dt = computeDt(nodes[neibNo], Rx, cellNo);
		if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
			traveltime =  nodes[neibNo].getTT(threadNo)+dt;
		}
	}
	return traveltime;
}


template<typename T1, typename T2>
	T1 Grid2Dvti_psv<T1,T2>::getTraveltime(const sxz<T1>& Rx,
																				const std::vector<Node2Dc<T1,T2> >& nodes,
																				T2& nodeParentRx, T2& cellParentRx,
																				const int threadNo) const {
    
	for ( size_t nn=0; nn<nodes.size(); ++nn ) {
		if ( nodes[nn] == Rx ) {
			nodeParentRx = nodes[nn].getNodeParent(threadNo);
			cellParentRx = nodes[nn].getCellParent(threadNo);
			return nodes[nn].getTT(threadNo);
		}
	}
    
	T2 cellNo = getCellNo( Rx );
	T2 neibNo = neighbors[cellNo][0];
	T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
	T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
	nodeParentRx = neibNo;
	cellParentRx = cellNo;
	for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
		neibNo = neighbors[cellNo][k];
		dt = computeDt(nodes[neibNo], Rx, cellNo);
		if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
			traveltime =  nodes[neibNo].getTT(threadNo)+dt;
			nodeParentRx = neibNo;
		}
	}
	return traveltime;
}







#endif
