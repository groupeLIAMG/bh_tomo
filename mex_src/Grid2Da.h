/*
 *  Grid2Da.h
 *  ttcr2da
 *
 *  Created by Bernard Giroux on 09-11-12.
 *  Copyright 2009 INRS ETE. All rights reserved.
 *
 */

#ifndef __GRID2DA_H__
#define __GRID2DA_H__

#include <stdint.h>

#include "Grid2Dc.h"

#ifdef _WIN32
#define pi 3.141592653589793
#define piD2 1.570796326794897
#endif

template<typename T1, typename T2>
class Grid2Da : public Grid2Dc<T1,T2> {
public:
	Grid2Da(const size_t nx, const size_t nz, const T1 ddx, const T1 ddz,
			const T1 minx, const T1 minz, const size_t nnx, const size_t nnz) :
	Grid2Dc<T1,T2>(nx, nz, ddx, ddz, minx, minz, nnx, nnz), xi(std::vector<T1>(nx*nz))
	{
    }

	T1 getSlowness(const size_t n) const {
		return 0.5 * Grid2Dc<T1,T2>::slowness[n] * ( 1.0 + std::sqrt(xi[n]) );
	}
	
	void setAnisotropyRatio(const T1 x) {
        for ( size_t n=0; n<xi.size(); ++n ) {
            xi[n] = x*x;
        }
    }
    
    int setAnisotropyRatio(const T1 *x, const size_t nx) {
        if ( xi.size() != nx ) {
            std::cerr << "Error: Anisotropy ratio vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<xi.size(); ++n ) {
            xi[n] = x[n]*x[n];
        }
        return 0;
    }
    
    int setAnisotropyRatio(const std::vector<T1>& x) {
        if ( xi.size() != x.size() ) {
            std::cerr << "Error: Anisotropy ratio vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<xi.size(); ++n ) {
            xi[n] = x[n]*x[n];
        }
        return 0;
    }
	
	
protected:
    std::vector<T1> xi;        // anisotropy ratio, xi = sz / sx, *** squared ***
    
	T1 computeDt(const Node2Dc<T1,T2>& source, const sxz<T1>& node,
				const T2 cellNo) const {
        T1 lx = node.x - source.getX();
        T1 lz = node.z - source.getZ();
        return Grid2Dc<T1,T2>::slowness[cellNo] *
		              std::sqrt( lx*lx + xi[cellNo]*lz*lz );
	}
	
	T1 computeDt(const Node2Dc<T1,T2>& source, const Node2Dc<T1,T2>& node,
				const T2 cellNo) const {
        T1 lx = node.getX() - source.getX();
        T1 lz = node.getZ() - source.getZ();
        return Grid2Dc<T1,T2>::slowness[cellNo] *
					std::sqrt( lx*lx + xi[cellNo]*lz*lz );
	}
	
	T1 get_tt_corr(const siv2<T1>& cell,
				  const Grid2Dc<T1,T2> *grid,
				  const size_t i) const {
		T1 t1 = Grid2Dc<T1,T2>::slowness[cell.i] * sqrt( cell.v1*cell.v1 +
					cell.v2*cell.v2 * xi[cell.i] );
		T1 t2 = sqrt( cell.v1*cell.v1 + cell.v2*cell.v2 ) * grid->getSlowness(i);
		return t1 - t2;
	}
	
	T1 getSlownessMoyenne(const subgridPar<T1,T2>& sgp) const {
		T1 smoy = 0.0;
		for ( size_t n=0; n<sgp.cells.size(); ++n ) {
			smoy += Grid2Dc<T1,T2>::slowness[ sgp.cells[n] ] *
			( 1.0 * xi[ sgp.cells[n] ] );
		}
		smoy /= 2*sgp.cells.size();
		return smoy;
	}

private:
	Grid2Da() {}
    Grid2Da(const Grid2Da<T1,T2>& g) {}
    Grid2Da<T1,T2>& operator=(const Grid2Da<T1,T2>& g) {}	
};




template<typename T1, typename T2>
class Grid2Daa : public Grid2Da<T1,T2> {
public:
	Grid2Daa(const size_t nx, const size_t nz, const T1 ddx, const T1 ddz,
			const T1 minx, const T1 minz, const size_t nnx, const size_t nnz) :
	Grid2Da<T1,T2>(nx, nz, ddx, ddz, minx, minz, nnx, nnz),
	aAngle(std::vector<T1>(nx*nz,0.0)), ca(std::vector<T1>(nx*nz,0.0)),
	sa(std::vector<T1>(nx*nz,0.0))
	{}
    
    T1 getAnisotropyAngle(const size_t n) const { return aAngle[n]; }
    
    void setAnisotropyAngle(const T1 s) {
        for ( size_t n=0; n<aAngle.size(); ++n ) {
            aAngle[n] = s;
			ca[n] = std::cos(s);
			sa[n] = std::sin(s);
        }
    }
    
    int setAnisotropyAngle(const T1 *s, const size_t ns) {
        if ( aAngle.size() != ns ) {
            std::cerr << "Error: anisotropy angle vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<aAngle.size(); ++n ) {
            aAngle[n] = s[n];
			ca[n] = std::cos(s[n]);
			sa[n] = std::sin(s[n]);
        }
        return 0;
    }
    
    int setAnisotropyAngle(const std::vector<T1>& s) {
        if ( aAngle.size() != s.size() ) {
            std::cerr << "Error: anisotropy angle vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<aAngle.size(); ++n ) {
            aAngle[n] = s[n];
			ca[n] = std::cos(s[n]);
			sa[n] = std::sin(s[n]);
        }
        return 0;
    }
    
protected:
#ifndef _WIN32
	static constexpr double pi = 3.141592653589793;
	static constexpr double piD2 = 1.570796326794897;
#endif

	std::vector<T1> aAngle;   // column-wise (z axis) anisotropy angle of the cells, in radians
	std::vector<T1> ca;       // cosine of aAngle
	std::vector<T1> sa;       // sine of aAngle
	
/*	double fast_atan2(const double y, const double x) const {
        
		if ( x == 0.0 ){
			if ( y > 0.0 ) return piD2;
			if ( y == 0.0 ) return 0.0;
			return -piD2;
		}
		
		double atan;
		double z = y/x;
		if ( fabs( z ) < 1.0 ){
			atan = z/(1.0 + 0.28*z*z);
			if ( x < 0.0 ){
				if ( y < 0.0 ) return atan - pi;
				return atan + pi;
			}
		}
		else{
			atan = piD2 - z/(z*z + 0.28f);
			if ( y < 0.0f ) return atan - pi;
		}
		return atan;
	}
	
	double fast_sin(const double x) const {
		double B = 4/pi;
		double C = -4/(pi*pi);
		double P = 0.225;
		
		double y = B * x + C * x * fabs(x);  //fast, inprecise
		
		// #ifdef EXTRA_PRECISION
		//  const float Q = 0.775;
		
		
		return y = P * (y * fabs(y) - y) + y;   // Q * y + P * y * fabs(y) //slower, precise
		// #endif
	}
	
	double fast_cos(const double x) const {
		double B = 4/pi;
		double C = -4/(pi*pi);
		double P = 0.225;
		double y;
		
		x = x + pi/2;
		if(x > pi){   // Original x > pi/2
			x -= 2 * pi;   // Wrap: cos(x) = cos(x - 2 pi)
		}
		
		y = B * x + C * x * fabs(x);  //fast, inprecise
		
		// #ifdef EXTRA_PRECISION
		//  const float Q = 0.775;
		
		
		return y = P * (y * fabs(y) - y) + y;   // Q * y + P * y * fabs(y) //slower, precise
		// #endif
	}*/
	
	
	T1 computeDt(const Node2Dc<T1,T2>& source, const sxz<T1>& node,
				const T2 cellNo) const {
		if ( aAngle[cellNo] == 0.0 ) {
            return Grid2Da<T1,T2>::computeDt(source, node, cellNo);
		}
		T1 lx = node.x - source.getX();
        T1 lz = node.z - source.getZ();

		T1 t1 = lx * ca[cellNo] + lz * sa[cellNo];
		T1 t2 = lz * ca[cellNo] - lx * sa[cellNo];
		
		return Grid2Dc<T1,T2>::slowness[cellNo] * 
		std::sqrt( t1*t1 + Grid2Da<T1,T2>::xi[cellNo]*Grid2Da<T1,T2>::xi[cellNo]*t2*t2 );
	}
	
	T1 computeDt(const Node2Dc<T1,T2>& source, const Node2Dc<T1,T2>& node,
				const T2 cellNo) const {
		if ( aAngle[cellNo] == 0.0 ) {
            return Grid2Da<T1,T2>::computeDt(source, node, cellNo);
		}
		T1 lx = node.getX() - source.getX();
        T1 lz = node.getZ() - source.getZ();
		
		T1 t1 = lx * ca[cellNo] + lz * sa[cellNo];
		T1 t2 = lz * ca[cellNo] - lx * sa[cellNo];
		
		return Grid2Dc<T1,T2>::slowness[cellNo] * 
		std::sqrt( t1*t1 + Grid2Da<T1,T2>::xi[cellNo]*Grid2Da<T1,T2>::xi[cellNo]*t2*t2 );
	}
	
private:
	Grid2Daa() {}
    Grid2Daa(const Grid2Daa<T1,T2>& g) {}
    Grid2Daa<T1,T2>& operator=(const Grid2Daa<T1,T2>& g) {}	
};

#endif
