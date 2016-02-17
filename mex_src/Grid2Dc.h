/*
 *  Grid2Dc.h
 *  ttRaisCourbes
 *
 *  Created by Bernard Giroux on 08-04-24.
 *  Copyright 2008 Bernard Giroux.
 *
 *
 * Reference paper
 *
 * @article{gruber:1062,
 *  author = {Thomas Gruber and Stewart A. Greenhalgh},
 *  collaboration = {},
 *  title = {Precision analysis of first-break times in grid models},
 *  publisher = {SEG},
 *  year = {1998},
 *  journal = {Geophysics},
 *  volume = {63},
 *  number = {3},
 *  pages = {1062-1065},
 *  url = {http://link.aip.org/link/?GPY/63/1062/1},
 *  doi = {10.1190/1.1444384}
 * }
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

#ifndef __GRID2DC_H__
#define __GRID2DC_H__




#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>

#include "Node2Dc.h"
#include "AntennaCorrection.h"




template<typename T1, typename T2>
struct subgridPar {
    subgridPar() : cells(std::vector<T2>()), nx(0), nz(0),
    ix_min(0), iz_min(0), dx(0), dz(0), xmin(0), zmin(0), nnx(0), nnz(0) {}
    std::vector<T2> cells;
    sxz<T1> c[5];
    T2 nx;
    T2 nz;
    T2 ix_min;
    T2 iz_min;
    T1 dx;
    T1 dz;
    T1 xmin;
    T1 zmin;
    T2 nnx;
    T2 nnz;
};


template<typename T1, typename T2>
class Grid2Dc {
public:
    Grid2Dc(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
           const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
		   const int nt=1);
    
    virtual ~Grid2Dc() {
        delete corr;
        delete subgridCellParent;
    }
    
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
    
    T2 getNumberOfCells() const { return slowness.size(); }
    virtual T1 getSlowness(const size_t n) const { return slowness[n]; }
    
    virtual void setSlowness(const T1 s) {
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s;
        }
    }
    
    virtual int setSlowness(const T1 *s, const size_t ns) {
        if ( slowness.size() != ns ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
        }
        return 0;
    }
    
    virtual int setSlowness(const std::vector<T1>& s) {
        if ( slowness.size() != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
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
                 std::vector<std::vector<siv<double> > >& l_data,
				 const int threadNo=0) const;
    
    int raytrace(const std::vector<sxz<T1> >& Tx,
                 const std::vector<T1>& t0, 
                 const std::vector<sxz<T1> >& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<std::vector<sxz<double> > >& r_data,
                 std::vector<std::vector<siv2<double> > >& l_data,
				 const int threadNo=0) const;
    
    int raytrace(const txPar<T1>& Tx,
                 const rxPar<T1>& Rx,
                 std::vector<T1>& traveltimes,
				 const int threadNo=0) const;
    
    int raytrace(const txPar<T1>& Tx,
                 const rxPar<T1>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<T1>& tt_corr,
                 std::vector<std::vector<sxz<double> > >& r_data,
                 std::vector<std::vector<siv<double> > >& l_data,
				 const int threadNo=0) const;
    
    int raytrace(const txPar<T1>& Tx,
                 const rxPar<T1>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<T1>& tt_corr,
                 std::vector<std::vector<sxz<double> > >& r_data,
                 std::vector<std::vector<siv2<double> > >& l_data,
				 const int threadNo=0) const;
    
    int buildAntCorr(const char *antenna ) {
        
        if ( (nsnx+1)%nsgx != 0 || (nsnz+1)%nsgz != 0 ) {
            std::cerr << "Error: sub grid number incorrect.\n";
            return 1;
        }
        // TODO: add other antennas!!!
        if ( corr != 0 ) delete corr;
        
        if ( strcmp(antenna, "Fixed 05")==0 )
            corr = new Fixed_05();
        else if ( strcmp(antenna, "Fixed 06")==0 )
            corr = new Fixed_06();
        else if ( strcmp(antenna, "Fixed 07")==0 )
            corr = new Fixed_07();
        else if ( strcmp(antenna, "Fixed 08")==0 )
            corr = new Fixed_08();
        else if ( strcmp(antenna, "Fixed 09")==0 )
            corr = new Fixed_09();
        else if ( strcmp(antenna, "Fixed 10")==0 )
            corr = new Fixed_10();
        else if ( strcmp(antenna, "Fixed 11")==0 )
            corr = new Fixed_11();
        else if ( strcmp(antenna, "Fixed 12")==0 )
            corr = new Fixed_12();
        else if ( strcmp(antenna, "Fixed 13")==0 )
            corr = new Fixed_13();
        else if ( strcmp(antenna, "Fixed 14")==0 )
            corr = new Fixed_14();
        else if ( strcmp(antenna, "Fixed 15")==0 )
            corr = new Fixed_15();
        else if ( strcmp(antenna, "Ramac 250")==0 )
            corr = new Ramac_250();
        else
            std::cout << "Warning: antenna correction undefined.  No correction applied\n";

        return 0;
    }
    
	
	int buildAntCorr(const T2 ngx, const T2 ngz, const char *antenna ) {
        nsgx = ngx;
        nsgz = ngz;

		return this->buildAntCorr( antenna );
	}
	
    void saveSlownessXYZ(const char filename[]) const {
        std::ofstream fout( filename );
        
        for ( T2 j=0, n=0; j<nCellx; ++j ) {
            T1 x = xmin + (0.5+j)*dx;
            for ( T2 i=0; i<nCellz; ++i, ++n ) {
                T1 z = zmin + (0.5+i)*dz;
                fout << x << "   " << z << "   " << slowness[n] << '\n';
            }
        }
        
        fout.close();
    }
    
    T2 getCellNo(const sxz<T1>& pt) const {
        T2 nx = static_cast<T2>( small + (pt.x-xmin)/dx );
        T2 nz = static_cast<T2>( small + (pt.z-zmin)/dz );
        return nx*nCellz + nz;
    }
    
    T2 getSubgridCellParent(const T2 i) const {
        return (*subgridCellParent)[i];
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
    
    std::vector<T1> slowness;   // column-wise (z axis) slowness vector of the cells
    std::vector<std::vector<T2> > neighbors;  // nodes common to a cell
    
    AntennaCorrection *corr;
    std::vector<T2> *subgridCellParent;
    
    void buildGridNodes();
    void buildGridNeighbors();
	
	AntennaCorrection* getAntennaCorrection() {return corr; }
    
    void propagate(std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
                   CompareNodePtr<T1> >& queue,
                   std::vector<bool>& inQueue,
                   std::vector<bool>& frozen,
				   const int threadNo) const;
	
	virtual T1 computeDt(const Node2Dc<T1,T2>& source, const sxz<T1>& node,
						const T2 cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
	}
    
	virtual T1 computeDt(const Node2Dc<T1,T2>& source, const Node2Dc<T1,T2>& node,
						const T2 cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
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
    
    void getSubGridPar(subgridPar<T1,T2>&, const sxz<T1>& pt, const T1 diam,
                       const T1 theta) const;
    Grid2Dc<T1,T2>* getSubGrid(const subgridPar<T1,T2>&, const T1 diam, const bool inWater) const;
    
    bool inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const;
    
    void attachSubgrid(Grid2Dc<T1,T2> *subGrid,
                       std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
                       CompareNodePtr<T1> >& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       bool toSubgrid) const;
    
    T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Dc<T1,T2> >& nodes,
					const int threadNo) const;
    
    T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Dc<T1,T2> >& nodes,
                    T2& nodeParentRx, T2& cellParentRx,
					const int threadNo) const;

    void save(const char filename[]) const;
    
	T1 get_tt_corr(const siv<T1>& cell,
				  const Grid2Dc<T1,T2> *grid,
				  const size_t i) const {
		return cell.v*(slowness[cell.i] - grid->slowness[i]);
	}
	
//	virtual T get_tt_corr(const siv2<T1>& cell,
//				  const Grid2Dc<T1,T2> *grid,
//				  const size_t i) const { // should never be used
//	}
	
	virtual T1 getSlownessMoyenne(const subgridPar<T1,T2>& sgp) const {
		T1 smoy = 0.0;
		for ( size_t n=0; n<sgp.cells.size(); ++n )
			smoy += slowness[ sgp.cells[n] ];
		smoy /= sgp.cells.size();
		return smoy;
	}
	
private:
    Grid2Dc() {}
    Grid2Dc(const Grid2Dc<T1,T2>& g) {}
    Grid2Dc<T1,T2>& operator=(const Grid2Dc<T1,T2>& g) {}
    
};

template<typename T1, typename T2>
Grid2Dc<T1,T2>::Grid2Dc(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                  const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
				  const int nt) : nThreads(nt),
dx(ddx), dz(ddz), xmin(minx), zmin(minz), xmax(minx+nx*ddx), zmax(minz+nz*ddz),
small(1.0e-10), nCellx(nx), nCellz(nz), nsnx(nnx), nsnz(nnz), nsgx(0), nsgz(0),
nodes(std::vector<Node2Dc<T1,T2> >( // noeuds secondaires
                              nCellx*nsnx*(nCellz+1) + nCellz*nsnz*(nCellx+1) +
                              // noeuds primaires
                              (nCellx+1) * (nCellz+1), Node2Dc<T1,T2>(nt) )),
slowness(std::vector<T1>(nCellx*nCellz)),
neighbors(std::vector<std::vector<T2> >(nCellx*nCellz)),
corr(0), subgridCellParent(0)
{
    buildGridNodes();
    buildGridNeighbors();
}

template<typename T1, typename T2>
void Grid2Dc<T1,T2>::buildGridNodes() {
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
            
            // primary nodes
            //            std::cout << n << "\t p\t-\t" << x << '\t' << z
            //            << "\t-\t" << cell_upLeft
            //            << '\t' << cell_downLeft
            //            << '\t' << cell_upRight
            //            << '\t' << cell_downRight << '\n';
            
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
                    
                    //                    std::cout << n << "\tsv\t-\t" << x << '\t' << zsv << "\t-\t" 
                    //                    << cell_downLeft << '\t' << cell_downRight << '\n';
                    
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
                    
                    //                    std::cout << n << "\tsh\t-\t" << xsh << '\t' << z << "\t-\t"
                    //                    << cell_upRight << '\t' << cell_downRight << '\n';
                    
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
            //            std::cout << '\n';
        }
        //        std::cout << '\n';
    }
}

template<typename T1, typename T2>
void Grid2Dc<T1,T2>::buildGridNeighbors() {
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
        }
    }
}


template<typename T1, typename T2>
int Grid2Dc<T1,T2>::check_pts(const std::vector<sxz<T1> >& pts) const {
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
void Grid2Dc<T1,T2>::initQueue(const std::vector<sxz<T1> >& Tx,
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
int Grid2Dc<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
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
int Grid2Dc<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
                        const std::vector<T1>& t0, 
                        const std::vector<sxz<T1> >& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxz<double> > >& r_data,
                        std::vector<std::vector<siv<double> > >& l_data,
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
    if ( l_data.size() != Rx.size() ) {
        l_data.resize( Rx.size() );
    }
    for ( size_t ni=0; ni<l_data.size(); ++ni ) {
        l_data[ni].resize( 0 );
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
        siv<double> cell;
		
        // store the son's coord 
        child.x = Rx[n].x;
        child.z = Rx[n].z;
        cell.i = cellParentRx;
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
 			
			r_tmp.push_back( child );

			cell.v = (*node_p)[iParent].getDistance( child );
			bool found=false;
			for (size_t nc=0; nc<l_data[n].size(); ++nc) {
				if ( l_data[n][nc].i == cell.i ) {
					l_data[n][nc].v += cell.v;
					found = true;
				}
			}
			if ( found == false ) {
				l_data[n].push_back( cell );
			}
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            cell.i = (*node_p)[iChild].getCellParent(threadNo);
			
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
		
		cell.v = (*node_p)[iParent].getDistance( child );
		bool found=false;
		for (size_t nc=0; nc<l_data[n].size(); ++nc) {
			if ( l_data[n][nc].i == cell.i ) {
				l_data[n][nc].v += cell.v;
				found = true;
			}
		}
		if ( found == false ) {
			l_data[n].push_back( cell );
		}
		
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );
		
        //  must be sorted to build matrix L
        sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
        
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
int Grid2Dc<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
                        const std::vector<T1>& t0, 
                        const std::vector<sxz<T1> >& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxz<double> > >& r_data,
                        std::vector<std::vector<siv2<double> > >& l_data,
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
    if ( l_data.size() != Rx.size() ) {
        l_data.resize( Rx.size() );
    }
    for ( size_t ni=0; ni<l_data.size(); ++ni ) {
        l_data[ni].resize( 0 );
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
        siv2<double> cell;
		
        // store the son's coord
        child.x = Rx[n].x;
        child.z = Rx[n].z;
        cell.i = cellParentRx;
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
			
			r_tmp.push_back( child );
			
			cell.v1 = (*node_p)[iParent].getDistanceX( child );
            cell.v2 = (*node_p)[iParent].getDistanceZ( child );
			bool found=false;
			for (size_t nc=0; nc<l_data[n].size(); ++nc) {
				if ( l_data[n][nc].i == cell.i ) {
					l_data[n][nc].v1 += cell.v1;
					l_data[n][nc].v2 += cell.v2;
					found = true;
				}
			}
			if ( found == false ) {
				l_data[n].push_back( cell );
			}
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            cell.i = (*node_p)[iChild].getCellParent(threadNo);

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
		
		cell.v1 = (*node_p)[iParent].getDistanceX( child );
		cell.v2 = (*node_p)[iParent].getDistanceZ( child );
		bool found=false;
		for (size_t nc=0; nc<l_data[n].size(); ++nc) {
			if ( l_data[n][nc].i == cell.i ) {
				l_data[n][nc].v1 += cell.v1;
				l_data[n][nc].v2 += cell.v2;
				found = true;
			}
		}
		if ( found == false ) {
			l_data[n].push_back( cell );
		}
		
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );

		//  must be sorted to build matrix L
        sort(l_data[n].begin(), l_data[n].end(), CompareSiv2_i<T1>());
        
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
int Grid2Dc<T1,T2>::raytrace(const txPar<T1>& Tx,
                        const rxPar<T1>& Rx,
                        std::vector<T1>& traveltimes,
						const int threadNo) const {
    
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    subgridPar<T1,T2> sgp;
    getSubGridPar(sgp, Tx.pt, Tx.diam, Tx.theta);
    Grid2Dc<T1,T2> *txGrid = getSubGrid(sgp, Tx.diam, Tx.inWater);
    if ( txGrid == 0 ) return 1;
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    std::vector<Node2Dc<T1,T2> > txNodes;
    std::vector<bool> inQueueTx( txGrid->nodes.size(), false );
    std::vector<bool> frozenTx( txGrid->nodes.size(), false );
    
    std::vector<sxz<T1> > txPt(1);
    txPt[0] = Tx.pt;
    std::vector<T1> t0(1, Tx.t0);
    
    txGrid->initQueue(txPt, t0, queue, txNodes, inQueueTx, frozenTx, threadNo);
    txGrid->propagate(queue, inQueueTx, frozenTx, threadNo);
    
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    attachSubgrid(txGrid, queue, inQueue, frozen, false);
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.pts.size() ) {
        traveltimes.resize( Rx.pts.size() );
    }
    
    for ( size_t nr=0; nr<Rx.pts.size(); ++nr ) {
        getSubGridPar(sgp, Rx.pts[nr], Rx.diam[nr], Rx.theta[nr]);
        Grid2Dc<T1,T2> *rxGrid = getSubGrid(sgp, Rx.diam[nr], Rx.inWater[nr]);
        if ( rxGrid == 0 ) { delete txGrid;  return 1; }
        
        std::vector<bool> inQueueRx( rxGrid->nodes.size(), false );
        std::vector<bool> frozenRx( rxGrid->nodes.size(), false );
        attachSubgrid(rxGrid, queue, inQueueRx, frozenRx, true);
        rxGrid->propagate(queue, inQueueRx, frozenRx, threadNo);
        
        traveltimes[nr] = rxGrid->getTraveltime(Rx.pts[nr], rxGrid->nodes,
												threadNo);
        delete rxGrid;
    }
    
    delete txGrid;
    return 0;
}

template<typename T1, typename T2>
int Grid2Dc<T1,T2>::raytrace(const txPar<T1>& Tx,
                        const rxPar<T1>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<T1>& tt_corr,
                        std::vector<std::vector<sxz<double> > >& r_data,
                        std::vector<std::vector<siv<double> > >& l_data,
						const int threadNo) const {//,
//                        const char *basename) const {
    
//	this->save("Grid2D.dat");
	
//	std::ofstream fout( "txPar.dat" );
//	fout << Tx.pt.x << ' ' << Tx.pt.z << ' ' << Tx.t0 << ' ' << Tx.theta << ' '
//	<< Tx.diam << ' ' << Tx.inWater;
//	fout.close();
//	
//	fout.open( "rxPar.dat ");
//	fout << Rx.pts.size();
//	for (size_t n=0; n<Rx.pts.size(); ++n)
//		fout << ' ' << Rx.pts[n].x << ' ' << Rx.pts[n].z;
//	fout << '\n';
//	fout << Rx.theta.size();
//	for (size_t n=0; n<Rx.theta.size(); ++n)
//		fout << ' ' << Rx.theta[n];
//	fout << '\n';
//	fout << Rx.diam.size();
//	for (size_t n=0; n<Rx.diam.size(); ++n)
//		fout << ' ' << Rx.diam[n];
//	fout << '\n';
//	fout << Rx.inWater.size();
//	for (size_t n=0; n<Rx.inWater.size(); ++n)
//		fout << ' ' << Rx.inWater[n];
//	fout << '\n';
//	fout.close();
	
	
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    subgridPar<T1,T2> sgp;
    getSubGridPar(sgp, Tx.pt, Tx.diam, Tx.theta);
    Grid2Dc<T1,T2> *txGrid = getSubGrid(sgp, Tx.diam, Tx.inWater);
    if ( txGrid == 0 ) return 1;
    
    //txGrid->saveSlownessXYZ( "txSubgrid.xyz" );
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    std::vector<Node2Dc<T1,T2> > txNodes;
    std::vector<bool> inQueueTx( txGrid->nodes.size(), false );
    std::vector<bool> frozenTx( txGrid->nodes.size(), false );
    
    std::vector<sxz<T1> > txPt(1);
    txPt[0] = Tx.pt;
    std::vector<T1> t0(1, Tx.t0);
    
    txGrid->initQueue(txPt, t0, queue, txNodes, inQueueTx, frozenTx, threadNo);
    txGrid->propagate(queue, inQueueTx, frozenTx, threadNo);
    
    
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    attachSubgrid(txGrid, queue, inQueue, frozen, false);
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.pts.size() ) {
        traveltimes.resize( Rx.pts.size() );
    }
    if ( tt_corr.size() != Rx.pts.size() ) {
        tt_corr.resize( Rx.pts.size() );
    }
    if ( l_data.size() != Rx.pts.size() ) {
        l_data.resize( Rx.pts.size() );
    }
    for ( size_t ni=0; ni<l_data.size(); ++ni ) {
        l_data[ni].resize( 0 );
    }
    if ( r_data.size() != Rx.pts.size() ) {
        r_data.resize( Rx.pts.size() );
    }
    for ( size_t ni=0; ni<r_data.size(); ++ni ) {
        r_data[ni].resize( 0 );
    }
    T2 nodeParentRx;
    T2 cellParentRx;
    
    
    for ( size_t nr=0; nr<Rx.pts.size(); ++nr ) {
        
        getSubGridPar(sgp, Rx.pts[nr], Rx.diam[nr], Rx.theta[nr]);
        Grid2Dc<T1,T2> *rxGrid = getSubGrid(sgp, Rx.diam[nr], Rx.inWater[nr]);
        
        if ( rxGrid == 0 ) { delete txGrid; return 1; }
        //rxGrid->saveSlownessXYZ( "rxSubgrid.xyz" );
        
        std::vector<bool> inQueueRx( rxGrid->nodes.size(), false );
        std::vector<bool> frozenRx( rxGrid->nodes.size(), false );
        attachSubgrid(rxGrid, queue, inQueueRx, frozenRx, true);
        rxGrid->propagate(queue, inQueueRx, frozenRx, threadNo);
        
        traveltimes[nr] = rxGrid->getTraveltime(Rx.pts[nr], rxGrid->nodes,
                                                nodeParentRx, cellParentRx,
												threadNo);
        tt_corr[nr] = 0.0;
        
        // Rx are in rxGrid->nodes
        //
        //  Start with Rx subgrid
        //
        
        std::vector<Node2Dc<T1,T2> > *node_p;
        node_p = &(rxGrid->nodes);
        
        std::vector<sxz<double> > r_tmp;
        std::vector<siv<double> > l_tmp;
        T2 iChild, iParent = nodeParentRx;
        sxz<double> child;
        siv<double> cell;

        // store the son's coord
        child.x = Rx.pts[nr].x;
        child.z = Rx.pts[nr].z;
        cell.i = rxGrid->getSubgridCellParent( cellParentRx );
        T2 iCellSubgrid = cellParentRx;
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {

			r_tmp.push_back( child );

            cell.v = (*node_p)[iParent].getDistance( child );
            tt_corr[nr] += get_tt_corr(cell, rxGrid, iCellSubgrid);
			
			bool found=false;
			for (size_t nc=0; nc<l_tmp.size(); ++nc) {
				if ( l_tmp[nc].i == cell.i ) {
					l_tmp[nc].v += cell.v;
					found = true;
				}
			}
			if ( found == false ) {
				l_tmp.push_back( cell );
			}
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            iCellSubgrid = (*node_p)[iChild].getCellParent(threadNo);
            cell.i = rxGrid->getSubgridCellParent( iCellSubgrid );

			iParent = (*node_p)[iParent].getNodeParent(threadNo);
        }
		// child is one node before arriving at the main grid

		r_tmp.push_back( child );
		
		cell.v = (*node_p)[iParent].getDistance( child );
		tt_corr[nr] += get_tt_corr(cell, rxGrid, iCellSubgrid);
		
		bool found=false;
		for (size_t nc=0; nc<l_tmp.size(); ++nc) {
			if ( l_tmp[nc].i == cell.i ) {
				l_tmp[nc].v += cell.v;
				found = true;
			}
		}
		if ( found == false ) {
			l_tmp.push_back( cell );
		}
		
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		
        //
        // Go on with main grid
        //
        // we are now at child.x - child.z
        node_p = &(nodes);
        for ( size_t nn=0; nn<node_p->size(); ++nn ) {
            if ( (*node_p)[nn] == child ) {
                iParent = nn;
                break;
            }
        }
        cell.i = (*node_p)[iParent].getCellParent(threadNo);
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
			
            r_tmp.push_back( child );
			
            cell.v = (*node_p)[iParent].getDistance( child );
            bool found=false;
			for (T2 nc=0; nc<l_tmp.size(); ++nc) {
				if ( l_tmp[nc].i == cell.i ) {
					l_tmp[nc].v += cell.v;
					found = true;
				}
			}
			if ( found == false ) {
				l_tmp.push_back( cell );
			}
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
			cell.i = (*node_p)[iChild].getCellParent(threadNo);
			
			iParent = (*node_p)[iParent].getNodeParent(threadNo);
        }
        
        // child is one node before arriving at the Tx subgrid
        
        r_tmp.push_back( child );
        
        cell.v = (*node_p)[iParent].getDistance( child );
        found=false;
        for (size_t nc=0; nc<l_tmp.size(); ++nc) {
            if ( l_tmp[nc].i == cell.i ) {
                l_tmp[nc].v += cell.v;
                found = true;
            }
        }
        if ( found == false ) {
            l_tmp.push_back( cell );
        }
        
        // we now go up in time - parent becomes the child of grand'pa
        iChild = iParent;
        child.x = (*node_p)[iChild].getX();
        child.z = (*node_p)[iChild].getZ();
        cell.i = (*node_p)[iChild].getCellParent(threadNo);

        //
        // End with Tx subgrid
        //
        // now child is on the Tx subgrid
        node_p = &(txGrid->nodes);
        for ( size_t nn=0; nn<node_p->size(); ++nn ) {
            if ( (*node_p)[nn] == child ) {
                iParent = nn;
                break;
            }
        }
        iCellSubgrid = (*node_p)[iParent].getCellParent(threadNo);
        cell.i = txGrid->getSubgridCellParent( iCellSubgrid );
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
			
			r_tmp.push_back( child );
			
            cell.v = (*node_p)[iParent].getDistance( child );
            tt_corr[nr] += get_tt_corr(cell, txGrid, iCellSubgrid);

			bool found=false;
			for (size_t nc=0; nc<l_tmp.size(); ++nc) {
				if ( l_tmp[nc].i == cell.i ) {
					l_tmp[nc].v += cell.v;
					found = true;
				}
			}
			if ( found == false ) {
				l_tmp.push_back( cell );
			}
            
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            iCellSubgrid = (*node_p)[iChild].getCellParent(threadNo);
            cell.i = txGrid->getSubgridCellParent( iCellSubgrid );

            iParent = (*node_p)[iParent].getNodeParent(threadNo);

            if ( iParent >= txGrid->nodes.size() ) {
                node_p = &txNodes;
                iParent -= txGrid->nodes.size();
            }
            else {
                node_p = &(txGrid->nodes);
            }
        }
		// parent is now at Tx
		r_tmp.push_back( child );
		
		cell.v = (*node_p)[iParent].getDistance( child );
		tt_corr[nr] += get_tt_corr(cell, txGrid, iCellSubgrid);
		
		found=false;
		for (size_t nc=0; nc<l_tmp.size(); ++nc) {
			if ( l_tmp[nc].i == cell.i ) {
				l_tmp[nc].v += cell.v;
				found = true;
			}
		}
		if ( found == false ) {
			l_tmp.push_back( cell );
		}
		
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );
		
		
        // the order should be from Tx to Rx, so we reorder...
        iParent = r_tmp.size();
        r_data[nr].resize( r_tmp.size() );
        for ( size_t nn=0; nn<r_data[nr].size(); ++nn ) {
            r_data[nr][nn].x = r_tmp[ iParent-1-nn ].x;
            r_data[nr][nn].z = r_tmp[ iParent-1-nn ].z;
        }
        
        sort(l_tmp.begin(), l_tmp.end(), CompareSiv_i<T1>());
        
        l_data[nr].push_back( l_tmp[0] );
        for ( size_t n1=1; n1<l_tmp.size(); ++n1 ) {
            bool found = false;
            for ( size_t n2=0; n2<l_data[nr].size(); ++n2) {
                if ( l_tmp[n1].i == l_data[nr][n2].i ) {
                    found = true;
                    l_data[nr][n2].v += l_tmp[n1].v;
                }
            }
            if ( !found ) {
                l_data[nr].push_back( l_tmp[n1] );
            }
        }
        /*
        if ( basename != 0 ) {
            char filename[100], val[100];
            sprintf(filename, "%s_grid_%02d.dat", basename, nr);
            std::ofstream fout( filename );
            for ( size_t n=0; n<nodes.size(); ++n ) {
                sprintf(val, "%lf  %lf  %10.9e\n", nodes[n].getX(), nodes[n].getZ(),
                        nodes[n].getTT());
                fout << val;
            }
            fout.close();
            
            sprintf(filename, "%s_sgrid_Tx_%02d.dat", basename, nr);
            fout.open( filename );
            for ( size_t n=0; n<txGrid->nodes.size(); ++n ) {
                sprintf(val, "%lf  %lf  %10.9e\n", txGrid->nodes[n].getX(),
                        txGrid->nodes[n].getZ(), txGrid->nodes[n].getTT());
                fout << val;
            }
            fout.close();
            
            sprintf(filename, "%s_sgrid_Rx_%02d.dat", basename, nr);
            fout.open( filename );
            for ( size_t n=0; n<rxGrid->nodes.size(); ++n ) {
                sprintf(val, "%lf  %lf  %10.9e\n", rxGrid->nodes[n].getX(),
                        rxGrid->nodes[n].getZ(), rxGrid->nodes[n].getTT());
                fout << val;
            }
            fout.close();
        }
        */
        delete rxGrid;
    }
    
    delete txGrid;
    return 0;
}

template<typename T1, typename T2>
int Grid2Dc<T1,T2>::raytrace(const txPar<T1>& Tx,
                        const rxPar<T1>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<T1>& tt_corr,
                        std::vector<std::vector<sxz<double> > >& r_data,
                        std::vector<std::vector<siv2<double> > >& l_data,
						const int threadNo) const {	
	
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    subgridPar<T1,T2> sgp;
    getSubGridPar(sgp, Tx.pt, Tx.diam, Tx.theta);
    Grid2Dc<T1,T2> *txGrid = getSubGrid(sgp, Tx.diam, Tx.inWater);
    if ( txGrid == 0 ) return 1;
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    std::vector<Node2Dc<T1,T2> > txNodes;
    std::vector<bool> inQueueTx( txGrid->nodes.size(), false );
    std::vector<bool> frozenTx( txGrid->nodes.size(), false );
    
    std::vector<sxz<T1> > txPt(1);
    txPt[0] = Tx.pt;
    std::vector<T1> t0(1, Tx.t0);
    
    txGrid->initQueue(txPt, t0, queue, txNodes, inQueueTx, frozenTx, threadNo);
    txGrid->propagate(queue, inQueueTx, frozenTx, threadNo);
    
    
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    attachSubgrid(txGrid, queue, inQueue, frozen, false);
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.pts.size() ) {
        traveltimes.resize( Rx.pts.size() );
    }
    if ( tt_corr.size() != Rx.pts.size() ) {
        tt_corr.resize( Rx.pts.size() );
    }
    if ( l_data.size() != Rx.pts.size() ) {
        l_data.resize( Rx.pts.size() );
    }
    for ( size_t ni=0; ni<l_data.size(); ++ni ) {
        l_data[ni].resize( 0 );
    }
    if ( r_data.size() != Rx.pts.size() ) {
        r_data.resize( Rx.pts.size() );
    }
    for ( size_t ni=0; ni<r_data.size(); ++ni ) {
        r_data[ni].resize( 0 );
    }
    T2 nodeParentRx;
    T2 cellParentRx;
    
    for ( size_t nr=0; nr<Rx.pts.size(); ++nr ) {
        
        getSubGridPar(sgp, Rx.pts[nr], Rx.diam[nr], Rx.theta[nr]);
        Grid2Dc<T1,T2> *rxGrid = getSubGrid(sgp, Rx.diam[nr], Rx.inWater[nr]);
        
        if ( rxGrid == 0 ) { delete txGrid; return 1; }
        
        std::vector<bool> inQueueRx( rxGrid->nodes.size(), false );
        std::vector<bool> frozenRx( rxGrid->nodes.size(), false );
        attachSubgrid(rxGrid, queue, inQueueRx, frozenRx, true);
        rxGrid->propagate(queue, inQueueRx, frozenRx, threadNo);
        
        traveltimes[nr] = rxGrid->getTraveltime(Rx.pts[nr], rxGrid->nodes,
                                                nodeParentRx, cellParentRx,
												threadNo);
        tt_corr[nr] = 0.0;
        
        // Rx are in rxGrid->nodes
        //
        //  Start with Rx subgrid
        //
        
        std::vector<Node2Dc<T1,T2> > *node_p;
        node_p = &(rxGrid->nodes);
        
        std::vector<sxz<double> > r_tmp;
        std::vector<siv2<double> > l_tmp;
        T2 iChild, iParent = nodeParentRx;
        sxz<double> child;
        siv2<double> cell;
        
        // store the son's coord
        child.x = Rx.pts[nr].x;
        child.z = Rx.pts[nr].z;
        cell.i = rxGrid->getSubgridCellParent( cellParentRx );
        T2 iCellSubgrid = cellParentRx;
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
            
			r_tmp.push_back( child );
            
            cell.v1 = (*node_p)[iParent].getDistanceX( child );
            cell.v2 = (*node_p)[iParent].getDistanceZ( child );
            tt_corr[nr] += get_tt_corr(cell, rxGrid, iCellSubgrid);
			
			bool found=false;
			for (size_t nc=0; nc<l_tmp.size(); ++nc) {
				if ( l_tmp[nc].i == cell.i ) {
					l_tmp[nc].v1 += cell.v1;
					l_tmp[nc].v2 += cell.v2;
					found = true;
				}
			}
			if ( found == false ) {
				l_tmp.push_back( cell );
			}
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            iCellSubgrid = (*node_p)[iChild].getCellParent(threadNo);
            cell.i = rxGrid->getSubgridCellParent( iCellSubgrid );
            
			iParent = (*node_p)[iParent].getNodeParent(threadNo);
        }
		// child is one node before arriving at the main grid
        
		r_tmp.push_back( child );
		
		cell.v1 = (*node_p)[iParent].getDistanceX( child );
		cell.v2 = (*node_p)[iParent].getDistanceZ( child );
		tt_corr[nr] += get_tt_corr(cell, rxGrid, iCellSubgrid);
		
		bool found=false;
		for (size_t nc=0; nc<l_tmp.size(); ++nc) {
			if ( l_tmp[nc].i == cell.i ) {
				l_tmp[nc].v1 += cell.v1;
				l_tmp[nc].v2 += cell.v2;
				found = true;
			}
		}
		if ( found == false ) {
			l_tmp.push_back( cell );
		}
		
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		
        //
        // Go on with main grid
        //
        // we are now at child.x - child.z
        node_p = &(nodes);
        for ( size_t nn=0; nn<node_p->size(); ++nn ) {
            if ( (*node_p)[nn] == child ) {
                iParent = nn;
                break;
            }
        }
        cell.i = (*node_p)[iParent].getCellParent(threadNo);
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
			
            r_tmp.push_back( child );
			
            cell.v1 = (*node_p)[iParent].getDistanceX( child );
            cell.v2 = (*node_p)[iParent].getDistanceZ( child );
            bool found=false;
			for (size_t nc=0; nc<l_tmp.size(); ++nc) {
				if ( l_tmp[nc].i == cell.i ) {
					l_tmp[nc].v1 += cell.v1;
					l_tmp[nc].v2 += cell.v2;
					found = true;
				}
			}
			if ( found == false ) {
				l_tmp.push_back( cell );
			}
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
			cell.i = (*node_p)[iChild].getCellParent(threadNo);
			
			iParent = (*node_p)[iParent].getNodeParent(threadNo);
        }
        
        // child is one node before arriving at the Tx subgrid
        
        r_tmp.push_back( child );
        
        cell.v1 = (*node_p)[iParent].getDistanceX( child );
        cell.v2 = (*node_p)[iParent].getDistanceZ( child );
        found=false;
        for (size_t nc=0; nc<l_tmp.size(); ++nc) {
            if ( l_tmp[nc].i == cell.i ) {
                l_tmp[nc].v1 += cell.v1;
                l_tmp[nc].v2 += cell.v2;
                found = true;
            }
        }
        if ( found == false ) {
            l_tmp.push_back( cell );
        }
        
        // we now go up in time - parent becomes the child of grand'pa
        iChild = iParent;
        child.x = (*node_p)[iChild].getX();
        child.z = (*node_p)[iChild].getZ();
        cell.i = (*node_p)[iChild].getCellParent(threadNo);
        
        //
        // End with Tx subgrid
        //
        // now child is on the Tx subgrid
        node_p = &(txGrid->nodes);
        for ( size_t nn=0; nn<node_p->size(); ++nn ) {
            if ( (*node_p)[nn] == child ) {
                iParent = nn;
                break;
            }
        }
        iCellSubgrid = (*node_p)[iParent].getCellParent(threadNo);
        cell.i = txGrid->getSubgridCellParent( iCellSubgrid );
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
			
			r_tmp.push_back( child );
			
            cell.v1 = (*node_p)[iParent].getDistanceX( child );
            cell.v2 = (*node_p)[iParent].getDistanceZ( child );
            tt_corr[nr] += get_tt_corr(cell, txGrid, iCellSubgrid);
            
			bool found=false;
			for (size_t nc=0; nc<l_tmp.size(); ++nc) {
				if ( l_tmp[nc].i == cell.i ) {
					l_tmp[nc].v1 += cell.v1;
					l_tmp[nc].v2 += cell.v2;
					found = true;
				}
			}
			if ( found == false ) {
				l_tmp.push_back( cell );
			}
            
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            iCellSubgrid = (*node_p)[iChild].getCellParent(threadNo);
            cell.i = txGrid->getSubgridCellParent( iCellSubgrid );
            
            iParent = (*node_p)[iParent].getNodeParent(threadNo);
            
            if ( iParent >= txGrid->nodes.size() ) {
                node_p = &txNodes;
                iParent -= txGrid->nodes.size();
            }
            else {
                node_p = &(txGrid->nodes);
            }
        }
		// parent is now at Tx
		r_tmp.push_back( child );
		
		cell.v1 = (*node_p)[iParent].getDistanceX( child );
		cell.v2 = (*node_p)[iParent].getDistanceZ( child );
		tt_corr[nr] += get_tt_corr(cell, txGrid, iCellSubgrid);
		
		found=false;
		for (size_t nc=0; nc<l_tmp.size(); ++nc) {
			if ( l_tmp[nc].i == cell.i ) {
				l_tmp[nc].v1 += cell.v1;
				l_tmp[nc].v2 += cell.v2;
				found = true;
			}
		}
		if ( found == false ) {
			l_tmp.push_back( cell );
		}
		
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );
		
		
        // the order should be from Tx to Rx, so we reorder...
        iParent = r_tmp.size();
        r_data[nr].resize( r_tmp.size() );
        for ( size_t nn=0; nn<r_data[nr].size(); ++nn ) {
            r_data[nr][nn].x = r_tmp[ iParent-1-nn ].x;
            r_data[nr][nn].z = r_tmp[ iParent-1-nn ].z;
        }
        
        sort(l_tmp.begin(), l_tmp.end(), CompareSiv2_i<T1>());
        
        l_data[nr].push_back( l_tmp[0] );
        for ( size_t n1=1; n1<l_tmp.size(); ++n1 ) {
            bool found = false;
            for ( size_t n2=0; n2<l_data[nr].size(); ++n2) {
                if ( l_tmp[n1].i == l_data[nr][n2].i ) {
                    found = true;
                    l_data[nr][n2].v1 += l_tmp[n1].v1;
                    l_data[nr][n2].v2 += l_tmp[n1].v2;
                }
            }
            if ( !found ) {
                l_data[nr].push_back( l_tmp[n1] );
            }
        }
        delete rxGrid;
    }
    
    delete txGrid;
    return 0;
}

template<typename T1, typename T2>
void Grid2Dc<T1,T2>::getSubGridPar(subgridPar<T1,T2>& sgp, const sxz<T1>& pt,
                              const T1 diam, const T1 theta) const {
    T1 length = 1.0;
    if ( corr != 0 ) length = corr->getLength();
    // find which cells are occupied by the antenna
    sgp.c[0].x =  0.501*diam;  sgp.c[0].z =  0.501*length;
    sgp.c[1].x =  0.501*diam;  sgp.c[1].z = -0.501*length;
    sgp.c[2].x = -0.501*diam;  sgp.c[2].z = -0.501*length;
    sgp.c[3].x = -0.501*diam;  sgp.c[3].z =  0.501*length;
    sgp.c[4].x =  0.0;       sgp.c[4].z =  0.0;
    
    T1 cosT = cos(theta);
    T1 sinT = sin(theta);
    
    for ( size_t n=0; n<4; ++n ) {
        T1 tmp   = sgp.c[n].x*cosT + sgp.c[n].z*sinT;
        sgp.c[n].z  = sgp.c[n].x*-sinT + sgp.c[n].z*cosT;
        sgp.c[n].x  = tmp;
        
        sgp.c[n].x += pt.x;
        sgp.c[n].z += pt.z;
        
        if ( sgp.c[n].x < xmin || sgp.c[n].x > xmax ||
            sgp.c[n].z < zmin || sgp.c[n].z > zmax ) {
            std::cerr << "Error: rotated antenna point outside the grid.\n";
            std::cerr << "       Point: " << sgp.c[n].x << ", "  << sgp.c[n].z << '\n';
            std::cerr << "       X min max: " << xmin << ' ' << xmax << '\n';
            std::cerr << "       Z min max: " << zmin << ' ' << zmax << '\n';
            return;
        }
    }
    // closing polygon
    sgp.c[4].x = sgp.c[0].x;
    sgp.c[4].z = sgp.c[0].z;
    
    T2 ix_max, iz_max;
    sgp.ix_min = ix_max = static_cast<T2>( (sgp.c[0].x-xmin)/dx);
    sgp.iz_min = iz_max = static_cast<T2>( (sgp.c[0].z-zmin)/dz);
    
    for ( size_t n=1; n<4; ++n ) {
        T2 tmp = static_cast<T2>( (sgp.c[n].x-xmin)/dx);
        sgp.ix_min = sgp.ix_min < tmp ? sgp.ix_min : tmp;
        ix_max = ix_max > tmp ? ix_max : tmp;
        tmp = static_cast<T2>( (sgp.c[n].z-zmin)/dz);
        sgp.iz_min = sgp.iz_min < tmp ? sgp.iz_min : tmp;
        iz_max = iz_max > tmp ? iz_max : tmp;
    }
    
    sgp.cells.resize(0);
    for ( T2 i=sgp.ix_min; i<=ix_max; ++i ) {
        for ( T2 j=sgp.iz_min; j<=iz_max; ++j ) {
            sgp.cells.push_back( i*nCellz + j );
        }
    }
    
    sgp.nx = nsgx*(1+ix_max-sgp.ix_min);
    sgp.nz = nsgz*(1+iz_max-sgp.iz_min);
    sgp.dx = dx/nsgx;
    sgp.dz = dz/nsgz;
    sgp.xmin = xmin + sgp.ix_min*dx;
    sgp.zmin = zmin + sgp.iz_min*dz;
    sgp.nnx = (nsnx+1)/nsgx - 1;
    sgp.nnz = (nsnz+1)/nsgz - 1;
}

template<typename T1, typename T2>
Grid2Dc<T1,T2>* Grid2Dc<T1,T2>::getSubGrid(const subgridPar<T1,T2>& sgp, const T1 diam,
                                 const bool inWater) const {
    
    Grid2Dc<T1,T2>* grid = new Grid2Dc<T1,T2>(sgp.nx, sgp.nz, sgp.dx, sgp.dz,
                                    sgp.xmin, sgp.zmin, sgp.nnx, sgp.nnz);
    if ( grid == 0 )
        return 0;
    
    T1 smoy = getSlownessMoyenne( sgp );
	
    T1 kTx = ( smoy * 0.2998 )*( smoy * 0.2998 );
    T1 sTx = 0.0;
    if ( corr != 0 ) sTx = corr->getSlowness(kTx, diam, inWater);
    
    grid->subgridCellParent = new std::vector<T2>(sgp.nx*sgp.nz);
    sxz<T1> p;
    for ( T2 ix=0, is=0; ix<sgp.nx; ++ix ) {
        // indices of main grid cell
        T2 ix2 = sgp.ix_min + ix/nsgx;
        p.x = sgp.xmin + 0.5*sgp.dx + ix*sgp.dx;
        for ( T2 iz=0; iz<sgp.nz; ++iz, ++is ) {
            T2 iz2 = sgp.iz_min + iz/nsgz;
            (*(grid->subgridCellParent))[is] = ix2*nCellz+iz2;
            
            // center of subcell
            p.z = sgp.zmin + 0.5*sgp.dz + iz*sgp.dz;
            if ( inPolygon( p, sgp.c, 5 ) && corr != 0 ) {
                grid->slowness[is] = sTx;
            }
            else {
                grid->slowness[is] = getSlowness( ix2*nCellz+iz2 );
            }
        }
    }
    return grid;
}


template<typename T1, typename T2>
void Grid2Dc<T1,T2>::propagate( std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
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
bool Grid2Dc<T1,T2>::inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const {
    bool c = false;
    for (size_t i = 0, j = N-1; i < N; j = i++) {
        if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
             ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
            (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
            c = !c;
    }
    return c;
}

template<typename T1, typename T2>
void Grid2Dc<T1,T2>::attachSubgrid(Grid2Dc<T1,T2> *subgrid,
                              std::priority_queue<Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
                              CompareNodePtr<T1> >& queue,
                              std::vector<bool>& inQueue,
                              std::vector<bool>& frozen,
                              bool toSubgrid) const {
    
    std::vector<Node2Dc<T1,T2>*> sgridNodes;
    std::vector<Node2Dc<T1,T2>*> gridNodes;
    
    T1 minx = subgrid->getXmin();
    T1 maxx = subgrid->getXmax();
    T1 minz = subgrid->getZmin();
    T1 maxz = subgrid->getZmax();
    
//    std::ofstream fout("attach_sgrid.dat");
    for ( size_t n=0; n<subgrid->nodes.size(); ++n ) {
        if (((fabs(subgrid->nodes[n].getX()-minx) < small || 
              fabs(subgrid->nodes[n].getX()-maxx) < small) &&
             (subgrid->nodes[n].getZ()>=(minz-small) && 
              subgrid->nodes[n].getZ()<=(maxz+small))) ||
            ((fabs(subgrid->nodes[n].getZ()-minz) < small || 
              fabs(subgrid->nodes[n].getZ()-maxz) < small) &&
             (subgrid->nodes[n].getX()>(minx+small) && 
              subgrid->nodes[n].getX()<(maxx-small)))) {
            sgridNodes.push_back( &(subgrid->nodes[n]) );
//            fout << subgrid->nodes[n].getX() << ' ' <<
//            subgrid->nodes[n].getZ() << '\n';
        }
    }
//    fout.close();
//    fout.open("attach_grid.dat");
    for ( size_t n=0; n<nodes.size(); ++n ) {
        if (((fabs(nodes[n].getX()-minx) < small || 
             fabs(nodes[n].getX()-maxx) < small) &&
            (nodes[n].getZ()>=(minz-small) && 
             nodes[n].getZ()<=(maxz+small))) ||
            ((fabs(nodes[n].getZ()-minz) < small || 
              fabs(nodes[n].getZ()-maxz) < small) &&
             (nodes[n].getX()>(minx+small) && 
              nodes[n].getX()<(maxx-small)))) {
            gridNodes.push_back( &(nodes[n]) );
//            fout << nodes[n].getX() << ' ' << nodes[n].getZ() << '\n';
        }
    }
//    fout.close();
//    exit(0);
    
    for ( size_t n1=0; n1<sgridNodes.size(); ++n1 ) {
        for ( size_t n2=0; n2<gridNodes.size(); ++n2 ) {
            if (sgridNodes[n1]->getDistance( *(gridNodes[n2]) ) < small) {
                if ( toSubgrid ) {
                    sgridNodes[n1]->setTT( gridNodes[n2]->getTT() );
                    queue.push( sgridNodes[n1] );
                    inQueue[ sgridNodes[n1]->getGridIndex() ] = true;
                    frozen[ sgridNodes[n1]->getGridIndex() ] = true;
                }
                else {
                    gridNodes[n2]->setTT( sgridNodes[n1]->getTT() );
                    queue.push( gridNodes[n2] );
                    inQueue[ gridNodes[n2]->getGridIndex() ] = true;
                    frozen[ gridNodes[n2]->getGridIndex() ] = true;
                }
                break;
            }
        }
    }
}


template<typename T1, typename T2>
T1 Grid2Dc<T1,T2>::getTraveltime(const sxz<T1>& Rx,
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
T1 Grid2Dc<T1,T2>::getTraveltime(const sxz<T1>& Rx,
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
void Grid2Dc<T1,T2>::save(const char filename[]) const {
	std::ofstream fout( filename );
	
	fout << dx << ' ' << dz << ' ' << xmin << ' ' << zmin << ' ' << xmax << ' ' << zmax << '\n';
	fout << nCellx << ' ' << nCellz << ' ' << nsnx << ' ' << nsnz << ' ' << nsgx << ' ' << nsgz << '\n';
//	fout << nodes.size() << '\n';
//	for ( size_t n=0; n < nodes.size(); ++n ) {
//		fout << nodes[n].getTT() << ' ' << nodes[n].getX() << ' ' << nodes[n].getZ() << ' '
//		<< ' ' << nodes[n].getGridIndex() << ' ' << nodes[n].getNodeParent() << ' '
//		<< nodes[n].getCellParent();
//		for (size_t no=0; no < nodes[n].getOwners().size(); ++no ) {
//			fout << ' ' << nodes[n].getOwners()[no];
//		}
//		fout << '\n';
//	}
	fout << slowness.size() << '\n';
	for ( size_t n=0; n < slowness.size(); ++n ) {
		fout << slowness[n] << '\n';
	}
//	fout << neighbors.size() << '\n';
//	for ( size_t n=0; n < neighbors.size(); ++n ) {
//		fout << neighbors[n].size();
//		for ( size_t nn=0; nn < neighbors[n].size(); ++nn ) {
//			fout << ' ' << neighbors[n][nn];
//		}
//		fout << '\n';
//	}
	if ( corr != 0 ) {
		fout << "1 " << corr->getName() << "\n";
	}
	else {
		fout << "0\n";
	}
	if ( subgridCellParent != 0 ) {
		fout << subgridCellParent->size();
		for ( size_t n=0; n < subgridCellParent->size(); ++n ) {
			fout << ' ' << (*subgridCellParent)[n];
		}
		fout << '\n';
	}
	else {
		fout << "0\n";
	}
	
	fout.close();
}

#endif
