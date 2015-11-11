//
//  Grid3D.h
//  ttcr.v2
//
//  Created by Giroux Bernard on 12-08-15.
//  Copyright (c) 2012 INRS-ETE. All rights reserved.
//

#ifndef __GRID3D_H__
#define __GRID3D_H__

namespace Grid {
    size_t get_number_of_nodes3D(const size_t nx, const size_t ny,
                                 const size_t nz, const size_t nnx,
                                 const size_t nny, const size_t nnz) {
        return nx*nnx*((ny+1)*(nz+1)) + // // secondary nodes on the edges
        ny*nny*((nx+1)*(nz+1)) +
        nz*nnz*((nx+1)*(ny+1)) +
        // secondary nodes on the faces
        (nnx*nny)*(nx*ny*(nz+1))+
        (nnx*nnz)*(nx*nz*(ny+1))+
        (nny*nnz)*(ny*nz*(nx+1))+
        // primary nodes
        (nx+1) * (ny+1) * (nz+1);
    }
};

template<typename T1, typename T2>
class Grid3D {
public:
    Grid3D(const T2 nx, const T2 ny, const T2 nz,
           const T1 ddx, const T1 ddy, const T1 ddz,
           const T1 minx, const T1 miny, const T1 minz,
           const T2 nnx, const T2 nny, const T2 nnz,
           const int nt=1) :
    nThreads(nt),
    dx(ddx), dy(ddy), dz(ddz),
    xmin(minx), ymin(miny), zmin(minz),
    xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
    nCellx(nx), nCelly(ny), nCellz(nz),
    nsnx(nnx), nsny(nny), nsnz(nnz)
    { }
    
    virtual ~Grid3D() {}
    
    T1 getDx() const { return dx; }
    T1 getDy() const { return dy; }
    T1 getDz() const { return dz; }
    T1 getXmin() const { return xmin; }
    T1 getXmax() const { return xmax; }
    T1 getYmin() const { return ymin; }
    T1 getYmax() const { return ymax; }
    T1 getZmin() const { return zmin; }
    T1 getZmax() const { return zmax; }
    T2 getNcellx() const { return nCellx; }
    T2 getNcelly() const { return nCelly; }
    T2 getNcellz() const { return nCellz; }
    T2 getNsnx() const { return nsnx; }
    T2 getNsny() const { return nsny; }
    T2 getNsnz() const { return nsnz; }

    T2 getNumberOfCells() const { return nCellx*nCelly*nCellz; }
	
    
protected:
    int nThreads;		     // number of threads
    T1 dx;                   // cell size in x
    T1 dy;			         // cell size in y
    T1 dz;                   // cell size in z
    T1 xmin;                 // x origin of the grid
    T1 ymin;                 // y origin of the grid
    T1 zmin;                 // z origin of the grid
    T1 xmax;                 // x end of the grid
    T1 ymax;                 // y end of the grid
    T1 zmax;                 // z end of the grid
    T2 nCellx;               // number of cells in x
    T2 nCelly;               // number of cells in y
    T2 nCellz;               // number of cells in z
    T2 nsnx;                 // number of secondary nodes in x
    T2 nsny;                 // number of secondary nodes in y
    T2 nsnz;                 // number of secondary nodes in z

};


#endif
