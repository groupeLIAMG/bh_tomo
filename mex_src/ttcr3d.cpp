//
//  Travel times with curved rays, in 3D
//
//  

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

#include "Grid3Dc.h"

using namespace std;

struct gridPar {
  double xmin;
	double ymin;
  double zmin;
  double dx;
	double dy;
  double dz;
  size_t nx;
	size_t ny;
  size_t nz;
  size_t nsx;
	size_t nsy;
  size_t nsz;
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double        *slowness, *xmin, *ymin, *zmin;
	double        *dx, *dy, *dz, *nx_d, *ny_d, *nz_d;
	double        *nsx_d, *nsy_d, *nsz_d;
  double        *Tx, *Rx, *t_arr, *Lval, *rays_p;
  size_t        nSlowness, stmp, nTx, nRx;
  mwSize        number_of_dims, nLmax;
  const mwSize  *dim_array;
  mxArray       **Rays;
  mwIndex       *irL, *jcL;

  gridPar gp;

  /* Check for proper number of arguments */
    
  if (nrhs != 4) {
		mexErrMsgTxt("TTCR3D requires four input arguments.");
  } else if (nlhs > 3) {
		mexErrMsgTxt("TTCR3D has a maximum of three output argument.");
  }
    
  // ------------------------------------------------------
	//	 Slowness
	// ------------------------------------------------------
    
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
	
  // ------------------------------------------------------
	//	 grid structure
  // ------------------------------------------------------
  if(!mxIsStruct(prhs[1]))
		mexErrMsgTxt("Grid must be a structure.");
	
  xmin = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "xmin") ) );
  ymin = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "ymin") ) );
  zmin = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "zmin") ) );
  dx = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "dx") ) );
  dy = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "dy") ) );
  dz = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "dz") ) );
  nx_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nx") ) );
  ny_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "ny") ) );
  nz_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nz") ) );
  nsx_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nsx") ) );
  nsy_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nsy") ) );
  nsz_d = static_cast<double*>( mxGetPr( mxGetField(prhs[1], 0, "nsz") ) );
	
#ifdef _WIN32
  gp.nx = static_cast<size_t>( *nx_d + 1.e-5 );
  gp.ny = static_cast<size_t>( *ny_d + 1.e-5 );
  gp.nz = static_cast<size_t>( *nz_d + 1.e-5 );
  gp.nsx = static_cast<size_t>( *nsx_d + 1.e-5 );
  gp.nsy = static_cast<size_t>( *nsy_d + 1.e-5 );
  gp.nsz = static_cast<size_t>( *nsz_d + 1.e-5 );
#else
  gp.nx = static_cast<size_t>( lround( *nx_d ) );
  gp.ny = static_cast<size_t>( lround( *ny_d ) );
  gp.nz = static_cast<size_t>( lround( *nz_d ) );
  gp.nsx = static_cast<size_t>( lround( *nsx_d ) );
  gp.nsy = static_cast<size_t>( lround( *nsy_d ) );
  gp.nsz = static_cast<size_t>( lround( *nsz_d ) );
#endif
  
  gp.xmin = *xmin;
  gp.ymin = *ymin;
  gp.zmin = *zmin;
  gp.dx = *dx;
  gp.dy = *dy;
  gp.dz = *dz;
  
  if ( nSlowness != gp.nx*gp.ny*gp.nz ) {
		mexErrMsgTxt("Slowness size inconsistent with grid params.");
  }
    
  // ------------------------------------------------------
	//  Tx
	// ------------------------------------------------------
	
  if (!(mxIsDouble(prhs[2]))) {
		mexErrMsgTxt("Tx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[2]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Tx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[2]);
  nTx = static_cast<size_t>( dim_array[0] );
  if ( dim_array[1] != 3 ) {
		mexErrMsgTxt("Tx: matrix nTx by 3.");
  }
  Tx = static_cast<double*>( mxGetPr(prhs[2]) );
	
  // ------------------------------------------------------
	//	 Rx
  // ------------------------------------------------------
	
  if (!(mxIsDouble(prhs[3]))) {
		mexErrMsgTxt("Rx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[3]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Rx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[3]);
  nRx = static_cast<size_t>( dim_array[0] );
  if ( dim_array[1] != 3 ) {
		mexErrMsgTxt("Rx: matrix nRx by 3.");
  }
  Rx = static_cast<double*>( mxGetPr(prhs[3]) );

  if ( nTx != nRx ) {
		mexErrMsgTxt("nTx should be equal to nRx.");
  }

  // ------------------------------------------------------
	//	 Output variable
  // ------------------------------------------------------
    
  plhs[0] = mxCreateDoubleMatrix(nRx, 1, mxREAL);
  t_arr = mxGetPr(plhs[0]);
    
    
  // ------------------------------------------------------
	// Optional output variables
  // ------------------------------------------------------
    
  if ( nlhs >= 2 ) {
		// 2rd arg: rays.
		plhs[1] = mxCreateCellMatrix(nRx, 1);
		Rays = (mxArray **) mxCalloc(nRx, sizeof(mxArray *));
  }
  // 3nd arg: L matrix. created after calculations

  // ------------------------------------------------------
	// Creation of the raytracing grid
  // ------------------------------------------------------
    
	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	mexPrintf("Using %d threads\n", nThreads);
	omp_set_num_threads( nThreads );
#endif

  Grid3Dc<double,uint32_t> grid(gp.nx, gp.ny, gp.nz,
																gp.dx, gp.dy, gp.dz,
																gp.xmin, gp.ymin, gp.zmin,
																gp.nsx, gp.nsy, gp.nsz, nThreads);

  if ( grid.setSlowness( slowness, nSlowness ) == 1 ) {
		mexErrMsgTxt("Problem with slowness assignement.");
  }


	//
	//	Looking for redundants Tx pts
  //
  vector< vector<sxyz<double> > > vTx;
  vector< vector<size_t> > iTx;
  sxyz<double> sxyz_tmp;
  sxyz_tmp.x = Tx[0];
	sxyz_tmp.y = Tx[nTx];
  sxyz_tmp.z = Tx[2*nTx];
  vTx.push_back( vector<sxyz<double> >(1, sxyz_tmp) );
  iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to
  // current Tx
  for ( size_t ntx=1; ntx<nTx; ++ntx ) {
		sxyz_tmp.x = Tx[ntx];
		sxyz_tmp.y = Tx[ntx+nTx];
		sxyz_tmp.z = Tx[ntx+2*nTx];
		bool found = false;
        
		for ( size_t nv=0; nv<vTx.size(); ++nv ) {
			if ( vTx[nv][0].x==sxyz_tmp.x &&
					 vTx[nv][0].y==sxyz_tmp.y &&
					 vTx[nv][0].z==sxyz_tmp.z ) {
				found = true;
				iTx[nv].push_back( ntx ) ;
				break;
			}
		}
		if ( !found ) {
			vTx.push_back( vector<sxyz<double> >(1, sxyz_tmp) );
			iTx.push_back( vector<size_t>(1, ntx) );
		}
  }


  //
	// Looping over all non redundant Tx
  //
	
	vector<vector<double> > tt( vTx.size() );
  vector<vector<siv<double> > > L_data(nTx);
  vector<vector<vector<siv<double> > > > l_data( vTx.size() );
  vector<vector<vector<sxyz<double> > > > r_data( vTx.size() );

	int nv, ni;

  vector<double> t0(1, 0.0);
  vector<sxyz<double> > vRx;


#pragma omp parallel //for default(none) shared(grid, vTx, t0, r_data, L_data, nlhs, Rx, nRx, t_arr, iTx, tt, l_data)	private(nv, ni, vRx, sxyz_tmp)
  for ( nv=0; nv<vTx.size(); ++nv ) {
		
		vRx.resize( 0 );
		for ( ni=0; ni<iTx[nv].size(); ++ni ) {
			sxyz_tmp.x = Rx[ iTx[nv][ni] ];
			sxyz_tmp.y = Rx[ iTx[nv][ni]+nRx ];
			sxyz_tmp.z = Rx[ iTx[nv][ni]+2*nRx ];
			vRx.push_back( sxyz_tmp );
		}
		
		if ( nlhs >= 2 ) {
			if ( grid.raytrace(vTx[nv], t0, vRx, tt[nv], l_data[nv], r_data[nv]
												 , omp_get_thread_num()) == 1 ) {
				mexErrMsgTxt("Problem while raytracing.");
			}
		}
		else {
			if ( grid.raytrace(vTx[nv], t0, vRx, tt[nv],
												 omp_get_thread_num()) == 1 ) {
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
				Rays[ iTx[nv][ni] ] = mxCreateDoubleMatrix(npts, 3, mxREAL);
				rays_p = (double*) mxGetData(Rays[ iTx[nv][ni] ]);
				for ( size_t np=0; np<npts; ++np ) {
					rays_p[np] = r_data[nv][ni][np].x;
					rays_p[np+npts] = r_data[nv][ni][np].y;
					rays_p[np+2*npts] = r_data[nv][ni][np].z;
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

}
