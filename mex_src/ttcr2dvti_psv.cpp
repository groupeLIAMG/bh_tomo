/*
 * P-SV waves travel times with curved rays, in 2D VTI media
 *
 * ttcr2dvti_psv
 *
 *  Created by Bernard Giroux on 13-03-06.
 *  Copyright 2013 Bernard Giroux.
 *
 * Reference paper
 *
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

#include "Grid2Dvti.h"

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
  double        *Vp0, *Vs0, *epsilon, *delta;
	double        *xmin, *zmin, *dx, *dz, *nx_d, *nz_d, *nsx_d, *nsz_d;
  double        *Tx, *Rx, *t_arr, *rays_p;
  size_t        nCells, nCells2, stmp, nTx, nRx;
  mwSize        number_of_dims, nLmax;
  const mwSize  *dim_array;
  mxArray       **Rays;
	int phase;
    
  gridPar gp;

    
  /* Check for proper number of arguments */
    
  if (nrhs != 8) {
		mexErrMsgTxt("TTCR2DVTI_SH requires eight input arguments.");
  } else if (nlhs > 2) {
		mexErrMsgTxt("TTCR2DVTI_SH has a maximum of two output argument.");
  }
    
  /* ------------------------------------------------------
		 Vp0
		 ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[0]))) {
		mexErrMsgTxt("Vp0 must be double precision.");
  }
  Vp0 = static_cast<double*>( mxGetPr(prhs[0]) );
  number_of_dims = mxGetNumberOfDimensions(prhs[0]);
  if ( number_of_dims != 2 ) {
		mexErrMsgTxt("Vp0 must be a vector (nCells by 1).");
  }
  dim_array = mxGetDimensions(prhs[0]);
  nCells = static_cast<size_t>( dim_array[0] );
  stmp = static_cast<size_t>( dim_array[1] );
  if ( stmp != 1 ) {
		mexErrMsgTxt("Vp0 must be a vector (nCells by 1).");
  }
	
  /* ------------------------------------------------------
		 Vs0
		 ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[1]))) {
		mexErrMsgTxt("Vs0 must be double precision.");
  }
  Vs0 = static_cast<double*>( mxGetPr(prhs[1]) );
  number_of_dims = mxGetNumberOfDimensions(prhs[1]);
  if ( number_of_dims != 2 ) {
		mexErrMsgTxt("Vs0 must be a vector (nCells by 1).");
  }
  dim_array = mxGetDimensions(prhs[1]);
  nCells2 = static_cast<size_t>( dim_array[0] );
  stmp = static_cast<size_t>( dim_array[1] );
  if ( stmp != 1 ) {
		mexErrMsgTxt("Vs0 must be a vector (nCells by 1).");
  }
	if ( nCells != nCells2 ) {
		mexErrMsgTxt("Vs0 and Vp0 must be the same size.");
  }

  /* ------------------------------------------------------
		 epsilon
		 ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[2]))) {
		mexErrMsgTxt("delta must be double precision.");
  }
  epsilon = static_cast<double*>( mxGetPr(prhs[2]) );
  number_of_dims = mxGetNumberOfDimensions(prhs[2]);
  if ( number_of_dims != 2 ) {
		mexErrMsgTxt("epsilon must be a vector (nCells by 1).");
  }
  dim_array = mxGetDimensions(prhs[2]);
  nCells2 = static_cast<size_t>( dim_array[0] );
  stmp = static_cast<size_t>( dim_array[1] );
  if ( stmp != 1 ) {
		mexErrMsgTxt("epsilon must be a vector (nCells by 1).");
  }
	if ( nCells != nCells2 ) {
		mexErrMsgTxt("epsilon and Vp0 must be the same size.");
  }
	
  /* ------------------------------------------------------
		 delta
		 ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[3]))) {
		mexErrMsgTxt("delta must be double precision.");
  }
  delta = static_cast<double*>( mxGetPr(prhs[3]) );
  number_of_dims = mxGetNumberOfDimensions(prhs[3]);
  if ( number_of_dims != 2 ) {
		mexErrMsgTxt("delta must be a vector (nCells by 1).");
  }
  dim_array = mxGetDimensions(prhs[3]);
  nCells2 = static_cast<size_t>( dim_array[0] );
  stmp = static_cast<size_t>( dim_array[1] );
  if ( stmp != 1 ) {
		mexErrMsgTxt("delta must be a vector (nCells by 1).");
  }
	if ( nCells != nCells2 ) {
		mexErrMsgTxt("delta and Vp0 must be the same size.");
  }
	
  /* ------------------------------------------------------
		 grid structure
     ------------------------------------------------------ */
  if(!mxIsStruct(prhs[4]))
		mexErrMsgTxt("Grid must be a structure.");
	
  xmin = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "xmin") ) );
  zmin = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "zmin") ) );
  dx = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "dx") ) );
  dz = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "dz") ) );
  nx_d = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "nx") ) );
  nz_d = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "nz") ) );
  nsx_d = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "nsx") ) );
  nsz_d = static_cast<double*>( mxGetPr( mxGetField(prhs[4], 0, "nsz") ) );
	
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
  
 
  gp.xmin = *xmin;
  gp.zmin = *zmin;
  gp.dx = *dx;
  gp.dz = *dz;
  
  if ( nCells != gp.nx*gp.nz ) {
		mexErrMsgTxt("Vs0 size inconsistent with grid params.");
  }
    
  /* ------------------------------------------------------
		 Tx
		 ------------------------------------------------------ */
	
  if (!(mxIsDouble(prhs[5]))) {
		mexErrMsgTxt("Tx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[5]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Tx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[5]);
  nTx = static_cast<size_t>( dim_array[0] );
  if ( dim_array[1] != 2 ) {
		mexErrMsgTxt("Tx: matrix nTx by 2.");
  }
  Tx = static_cast<double*>( mxGetPr(prhs[5]) );
	
  /* ------------------------------------------------------
		 Rx
     ------------------------------------------------------ */
	
  if (!(mxIsDouble(prhs[6]))) {
		mexErrMsgTxt("Rx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[6]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Rx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[6]);
  nRx = static_cast<size_t>( dim_array[0] );
  if ( dim_array[1] != 2 ) {
		mexErrMsgTxt("Rx: matrix nRx by 2.");
  }
  Rx = static_cast<double*>( mxGetPr(prhs[6]) );

  if ( nTx != nRx ) {
		mexErrMsgTxt("nTx should be equal to nRx.");
  }

	/* ------------------------------------------------------
		 phase
     ------------------------------------------------------ */

  if (!(mxIsLogicalScalar(prhs[7]))) {
		mexErrMsgTxt("phase must be logical scalar.");
  }
	if ( mxIsLogicalScalarTrue(prhs[7]) ) {
		phase = 1;
	}
	else {
		phase = 0;
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

  /* ------------------------------------------------------
		 Creation of the raytracing grid
     ------------------------------------------------------ */
    
	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	mexPrintf("Using %d threads\n", nThreads);
	omp_set_num_threads( nThreads );
#endif
		
  Grid2Dvti_psv<double,uint32_t> grid(gp.nx, gp.nz, gp.dx, gp.dz, gp.xmin, gp.zmin,
																			gp.nsx, gp.nsz, nThreads);
    
  if ( grid.setVp0( Vp0, nCells ) == 1 ) {
		mexErrMsgTxt("Problem with Vp0 assignement.");
  }
  if ( grid.setVs0( Vs0, nCells ) == 1 ) {
		mexErrMsgTxt("Problem with Vs0 assignement.");
  }
	if ( grid.setEpsilon( epsilon, nCells ) == 1 ) {
		mexErrMsgTxt("Problem with epsilon assignement.");
  }
	if ( grid.setDelta( delta, nCells ) == 1 ) {
		mexErrMsgTxt("Problem with delta assignement.");
  }
	grid.setPhase( phase );
    
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
		
		if ( nlhs == 2 ) {
			if ( grid.raytrace(vTx[nv], t0, vRx, tt[nv], r_data[nv], omp_get_thread_num()) == 1 ) {
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

	if ( nlhs == 2 ) {
		for ( nv=0; nv<vTx.size(); ++nv ) {

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
	
 	
  //return;
}
