/*
 * Travel times with curved rays, in 2D
 *
 * ttcr2d 
 *
 *  Created by Bernard Giroux on 08-04-26.
 *  Copyright 2008 Bernard Giroux.
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
 */

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>

#ifdef _OPENMP
#include "omp.h"
#else
#define omp_get_thread_num() 0
#endif

#include "mex.h"

#include "Grid2Dc.h"

using namespace std;

extern void _main();

struct gridPar {
  double xmin;
  double zmin;
  double dx;
  double dz;
  size_t nx;
  size_t nz;
  size_t nsx;
  size_t nsz;
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double        *slowness, *xmin, *zmin, *dx, *dz, *nx_d, *nz_d, *nsx_d, *nsz_d;
  double        *Tx, *Rx, *t_arr, *Lval, *rays_p;
  size_t        nSlowness, stmp, nTx, nRx;
  mwSize        number_of_dims, nLmax;
  const mwSize  *dim_array;
  mxArray       **Rays;
  mwIndex       *irL, *jcL;
    
  gridPar gp;

    
  /* Check for proper number of arguments */
    
  if (nrhs != 4) {
		mexErrMsgTxt("TTCR2D requires four input arguments.");
  } else if (nlhs > 3) {
		mexErrMsgTxt("TTCR2D has a maximum of three output argument.");
  }
    
  /* ------------------------------------------------------
		 Slowness
		 ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[0]))) {
		mexErrMsgTxt("Slowness must be double precision.");
  }
  slowness = static_cast<double*>( mxGetPr(prhs[0]) );
  number_of_dims = mxGetNumberOfDimensions(prhs[0]);
  if ( number_of_dims != 2 ) {
		mexErrMsgTxt("Slowness must be a vector (nSlowness by 1).");
  }
  dim_array = mxGetDimensions(prhs[0]);
  nSlowness = static_cast<size_t>( dim_array[0] );
  stmp = static_cast<size_t>( dim_array[1] );
  if ( stmp != 1 ) {
		mexErrMsgTxt("Slowness must be a vector (nSlowness by 1).");
  }
	
  /* ------------------------------------------------------
		 grid structure
     ------------------------------------------------------ */    
  if(!mxIsStruct(prhs[1]))
		mexErrMsgTxt("Grid must be a structure.");
	
  xmin = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "xmin") ) );
  zmin = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "zmin") ) );
  dx = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "dx") ) );
  dz = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "dz") ) );
  nx_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nx") ) );
  nz_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nz") ) );
  nsx_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nsx") ) );
  nsz_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nsz") ) );
	
#ifdef _WIN32
  gp.nx = static_cast<size_t>( *nx_d + 1.e-5 );
  gp.nz = static_cast<size_t>( *nz_d + 1.e-5 );
  gp.nsx = static_cast<size_t>( *nsx_d + 1.e-5 );
  gp.nsz = static_cast<size_t>( *nsz_d + 1.e-5 );
#else
  gp.nx = static_cast<size_t>( lround( *nx_d ) );
  gp.nz = static_cast<size_t>( lround( *nz_d ) );
  gp.nsx = static_cast<size_t>( lround( *nsx_d ) );
  gp.nsz = static_cast<size_t>( lround( *nsz_d ) );
#endif
  
  //    cout << nSlowness << " - " << *xmin << " " << *zmin << " - "
  //    << *dx << " " << *dz << " - "
  //    << gp.nx << " " << gp.nz << " - " << gp.nsx << " " << gp.nsz << '\n';
  
  gp.xmin = *xmin;
  gp.zmin = *zmin;
  gp.dx = *dx;
  gp.dz = *dz;
  
  if ( nSlowness != gp.nx*gp.nz ) {
		mexErrMsgTxt("Slowness size inconsistent with grid params.");
  }
    
  /* ------------------------------------------------------
		 Tx
		 ------------------------------------------------------ */
	
  if (!(mxIsDouble(prhs[2]))) {
		mexErrMsgTxt("Tx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[2]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Tx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[2]);
  nTx = static_cast<size_t>( dim_array[0] );
  if ( dim_array[1] != 2 ) {
		mexErrMsgTxt("Tx: matrix nTx by 2.");
  }
  Tx = static_cast<double*>( mxGetPr(prhs[2]) );
	
  /* ------------------------------------------------------
		 Rx
     ------------------------------------------------------ */
	
  if (!(mxIsDouble(prhs[3]))) {
		mexErrMsgTxt("Rx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[3]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Rx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[3]);
  nRx = static_cast<size_t>( dim_array[0] );
  if ( dim_array[1] != 2 ) {
		mexErrMsgTxt("Rx: matrix nRx by 2.");
  }
  Rx = static_cast<double*>( mxGetPr(prhs[3]) );

  if ( nTx != nRx ) {
		mexErrMsgTxt("nTx should be equal to nRx.");
  }

  /* ------------------------------------------------------
		 Output variable
     ------------------------------------------------------ */
    
  plhs[0] = mxCreateDoubleMatrix(nRx, 1, mxREAL);
  t_arr = mxGetPr(plhs[0]);
    
    
  /* ------------------------------------------------------
		 Optional output variables
     ------------------------------------------------------ */
    
  if ( nlhs >= 2 ) {
		// 2rd arg: rays.
		plhs[1] = mxCreateCellMatrix(nRx, 1);
		Rays = (mxArray **) mxCalloc(nRx, sizeof(mxArray *));
  }
  // 3nd arg: L matrix. created after calculations

  /* ------------------------------------------------------
		 Creation of the raytracing grid
     ------------------------------------------------------ */
    
	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	mexPrintf("Using %d threads\n", nThreads);
	omp_set_num_threads( nThreads );
#endif
		
  Grid2Dc<double,uint32_t> grid(gp.nx, gp.nz, gp.dx, gp.dz, gp.xmin, gp.zmin,
																gp.nsx, gp.nsz, nThreads);
    
  if ( grid.setSlowness( slowness, nSlowness ) == 1 ) {
		mexErrMsgTxt("Problem with slowness assignement.");
  }
	
    
  /*
		Looking for redundants Tx pts
  */
  vector< vector<sxz<double> > > vTx;
  vector< vector<size_t> > iTx;
  sxz<double> sxz_tmp;
  sxz_tmp.x = Tx[0];
  sxz_tmp.z = Tx[nTx];
  vTx.push_back( vector<sxz<double> >(1, sxz_tmp) );
  iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to
  // current Tx
  for ( size_t ntx=1; ntx<nTx; ++ntx ) {
		sxz_tmp.x = Tx[ntx];
		sxz_tmp.z = Tx[ntx+nTx];
		bool found = false;
        
		for ( size_t nv=0; nv<vTx.size(); ++nv ) {
			if ( vTx[nv][0].x==sxz_tmp.x && vTx[nv][0].z==sxz_tmp.z ) {
				found = true;
				iTx[nv].push_back( ntx ) ;
				break;
			}
		}
		if ( !found ) {
			vTx.push_back( vector<sxz<double> >(1, sxz_tmp) );
			iTx.push_back( vector<size_t>(1, ntx) );
		}
  }
  vector<double> t0(1, 0.0);
	

  /*
		Looping over all non redundant Tx
  */
	
  vector<sxz<double> > vRx;
	vector<vector<double> > tt( vTx.size() );
  vector<vector<siv<double> > > L_data(nTx);
  vector<vector<vector<siv<double> > > > l_data( vTx.size() );
  vector<vector<vector<sxz<double> > > > r_data( vTx.size() );

	int nv, ni;
#pragma omp parallel for default(none) shared(grid, vTx, t0, r_data, L_data, nlhs, Rx, nRx, t_arr, iTx, tt, l_data)	private(nv, ni, vRx, sxz_tmp)
  for ( nv=0; nv<vTx.size(); ++nv ) {
		
		vRx.resize( 0 );
		for ( ni=0; ni<iTx[nv].size(); ++ni ) {
			sxz_tmp.x = Rx[ iTx[nv][ni] ];
			sxz_tmp.z = Rx[ iTx[nv][ni]+nRx ];
			vRx.push_back( sxz_tmp );
		}
		
		if ( nlhs >= 2 ) {
			if ( grid.raytrace(vTx[nv], t0, vRx, tt[nv], r_data[nv], l_data[nv], omp_get_thread_num()) == 1 ) {
				mexErrMsgTxt("Problem while raytracing.");
			}
		}
		else {
			if ( grid.raytrace(vTx[nv], t0, vRx, tt[nv], omp_get_thread_num()) == 1 ) {
				mexErrMsgTxt("Problem while raytracing.");
			}
		}
		
	} /*-- End of omp parallel for --*/

  for ( nv=0; nv<vTx.size(); ++nv ) {
		for ( ni=0; ni<iTx[nv].size(); ++ni ) {
			t_arr[ iTx[nv][ni] ] = tt[nv][ni];
		}
	}

	if ( nlhs >= 2 ) {
		for ( nv=0; nv<vTx.size(); ++nv ) {

			for ( ni=0; ni<iTx[nv].size(); ++ni ) {
				L_data[ iTx[nv][ni] ] = l_data[nv][ni];
			}

			for ( ni=0; ni<iTx[nv].size(); ++ni ) {
				size_t npts = r_data[nv][ni].size();
				Rays[ iTx[nv][ni] ] = mxCreateDoubleMatrix(npts, 2, mxREAL);
				rays_p = (double*) mxGetData(Rays[ iTx[nv][ni] ]);
				for ( size_t np=0; np<npts; ++np ) {
					rays_p[np] = r_data[nv][ni][np].x;
					rays_p[np+npts] = r_data[nv][ni][np].z;
				}
				mxSetCell( plhs[1], iTx[nv][ni], Rays[ iTx[nv][ni] ] );
			}
		}   
  }
	
  
  if ( nlhs == 3 ) {
		nLmax = 0;
		for ( size_t n=0; n<L_data.size(); ++n ) {
			nLmax += L_data[n].size();
		}
		plhs[2] = mxCreateSparse(nTx, nSlowness, nLmax, mxREAL);
		Lval = mxGetPr( plhs[2] );
		irL  = mxGetIr( plhs[2] );
		jcL  = mxGetJc( plhs[2] );
		
		size_t k = 0;
		for ( size_t j=0; j<nSlowness; ++j ) {
			jcL[j] = k;
			for ( size_t i=0; i<nTx; ++i ) {
				for ( size_t n=0; n<L_data[i].size(); ++n ) {
					if ( L_data[i][n].i == j ) {
						irL[k] = i;
						Lval[k] = L_data[i][n].v;
						k++;
					}
				}
			}
		}
		jcL[nSlowness] = k;
  }
	
  //return;
}
