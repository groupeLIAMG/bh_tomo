/*
 * Matrix L from straight rays, in 2D
 *
 * Lsr2d
 *
 *  Created by Bernard Giroux on 08-04-26.
 *  Copyright 2008 Bernard Giroux.
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

#include <math.h>

#include "mex.h"

extern void _main();



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double        dx, dz, xs, zs, xr, zr, dl, dlx, dlz, x, z, dtmp, m, b;
  double        percent_sp, xe, ze, zi;
  const double  small=1.e-10;
  double        *Tx, *Rx, *grx, *grz, *Lval, *gridx, *gridz, *Lvec;
  mwSize        nTx, nRx, n_grx, n_grz, nSlowness;
  mwSize        number_of_dims, nLmax;
  const mwSize  *dim_array;
  mwIndex       *irL, *jcL, *iVec, *jVec;
  mwSize        n, ix, iz, iCell, i, j, k, nL, oldnzmax;
	int           up;

    
  /* Check for proper number of arguments */
    
  if (nrhs != 4) {
		mexErrMsgTxt("LSR2D requires four input arguments.");
  } else if (nlhs > 3) {
		mexErrMsgTxt("LSR2D has a maximum of three output argument.");
  }
    
    
  /* ------------------------------------------------------
		 Tx
     ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[0]))) {
		mexErrMsgTxt("Tx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[0]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Tx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[0]);
  nTx = dim_array[0];
  if ( dim_array[1] != 2 ) {
		mexErrMsgTxt("Tx: matrix nTx by 2.");
  }
  Tx = (double*)mxGetPr(prhs[0]);
    
  /* ------------------------------------------------------
		 Rx
     ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[1]))) {
		mexErrMsgTxt("Rx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[1]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("Rx must be a rank 2 matrix.");
  }
  dim_array = mxGetDimensions(prhs[1]);
  nRx = dim_array[0];
  if ( dim_array[1] != 2 ) {
		mexErrMsgTxt("Rx: matrix nRx by 2.");
  }
  Rx = (double*)mxGetPr(prhs[1]);

  if ( nTx != nRx ) {
		mexErrMsgTxt("nTx should be equal to nRx.");
  }

  /* ------------------------------------------------------
		 grx
     ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[2]))) {
		mexErrMsgTxt("grx must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[2]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("grx must be a vector.");
  }
  dim_array = mxGetDimensions(prhs[2]);
  if ( dim_array[0] != 1 && dim_array[1] != 1 ) {
		mexErrMsgTxt("grx: vector 1 by n_grx.");
  }
  n_grx = dim_array[1] > dim_array[0] ? dim_array[1] : dim_array[0];
  grx = (double*)mxGetPr(prhs[2]);
  dx = grx[1]-grx[0];

  /* ------------------------------------------------------
		 grz
     ------------------------------------------------------ */
    
  if (!(mxIsDouble(prhs[3]))) {
		mexErrMsgTxt("grz must be double precision.");
  }
  number_of_dims = mxGetNumberOfDimensions(prhs[3]);
  if ( number_of_dims != 2 ){
		mexErrMsgTxt("grz must be a vector.");
  }
  dim_array = mxGetDimensions(prhs[3]);
  if ( dim_array[0] != 1 && dim_array[1] != 1 ) {
		mexErrMsgTxt("grz: vector 1 by n_grx.");
  }
  n_grz = dim_array[1] > dim_array[0] ? dim_array[1] : dim_array[0];
  grz = (double*)mxGetPr(prhs[3]);
  dz = grz[1]-grz[0];




  nSlowness = (n_grx-1)*(n_grz-1);
  nLmax = nTx * n_grx * n_grz/2;
  percent_sp = (nLmax*1.0)/(nTx*nSlowness*1.0);

  Lvec = mxMalloc( nLmax*sizeof(double) );
  iVec = mxMalloc( nLmax*sizeof(mwIndex) );
  jVec = mxMalloc( nLmax*sizeof(mwIndex) );

  k = 0;
  for ( n=0; n<nTx; ++n ) {

		xs = Tx[n];
		zs = Tx[nTx+n];
		xr = Rx[n];
		zr = Rx[nRx+n];

		if ( xs>xr ) {  /* on va de s à r, on veut x croissant */
			dtmp = xs;
			xs = xr;
			xr = dtmp;
			dtmp = zs;
			zs = zr;
			zr = dtmp;
		}

		/* points de depart */
		x = xs;
		z = zs;

		if ( fabs(zs-zr)<small ) {  /* rai horizontal */

			for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
			for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

			while ( x < xr ) {
				iCell = ix*(n_grz-1) + iz;

				dlx = ( grx[ix+1]<xr ? grx[ix+1] : xr ) - x;

				iVec[k] = n;
				jVec[k] = iCell;
				Lvec[k] = dlx;
				k++;

				if (k>=nLmax){
					oldnzmax = nLmax;
					percent_sp += 0.1;
					nLmax = (mwSize)ceil((double)nTx*(double)nSlowness*percent_sp);

					/* make sure nzmax increases at least by 1 */
					if (oldnzmax == nLmax) 
						nLmax++;
					mxRealloc( Lvec, nLmax*sizeof(double) );
					mxRealloc( iVec, nLmax*sizeof(mwIndex) );
					mxRealloc( jVec, nLmax*sizeof(mwIndex) );
				}

				ix++;
				x = grx[ix];
			}
		}
		else if ( fabs(xs-xr)<small ) { /* rai vertical */
			if ( zs > zr ) {  /* on va de s à r, on veut z croissant */
				dtmp = zs;
				zs = zr;
				zr = dtmp;
			}
			z = zs;

			for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
			for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

			while ( z < zr ) {
				iCell = ix*(n_grz-1) + iz;

				dlz = ( grz[iz+1]<zr ? grz[iz+1] : zr ) - z;

				iVec[k] = n;
				jVec[k] = iCell;
				Lvec[k] = dlz;
				k++;

				if (k>=nLmax){
					oldnzmax = nLmax;
					percent_sp += 0.1;
					nLmax = (mwSize)ceil((double)nTx*(double)nSlowness*percent_sp);

					/* make sure nzmax increases at least by 1 */
					if (oldnzmax == nLmax) 
						nLmax++;
					mxRealloc( Lvec, nLmax*sizeof(double) );
					mxRealloc( iVec, nLmax*sizeof(mwIndex) );
					mxRealloc( jVec, nLmax*sizeof(mwIndex) );
				}

				iz++;
				z = grz[iz];
			}
		}
		else { /* rai oblique */

			/* pente du rai */
			m = (zr-zs)/(xr-xs);
			b = zr - m*xr;
			up = m>0;

			for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
			for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

			while ( x < xr ) {
		
				zi = m*grx[ix+1] + b;

				if ( up ) {
					while ( z < zi && z < zr ) {
						iCell = ix*(n_grz-1) + iz;
			
						ze = grz[iz+1]<zi ? grz[iz+1] : zi;
						ze = ze<zr ? ze : zr;
						xe = (ze-b)/m;
						dlx = xe - x;
						dlz = ze - z;
						dl = sqrt( dlx*dlx + dlz*dlz );
			
						/* mexPrintf("x %lf xe %lf      z %lf ze %lf  zi %lf\n", x, xe, z, ze, zi); */

						/* mexPrintf("i = %d, j = %d, dl = %lf, dlx = %lf, dlz = %lf\n", ix, iz, dl, dlx, dlz); */
						iVec[k] = n;
						jVec[k] = iCell;
						Lvec[k] = dl;
						k++;
			
						if (k>=nLmax){
							oldnzmax = nLmax;
							percent_sp += 0.1;
							nLmax = (mwSize)ceil((double)nTx*(double)nSlowness*percent_sp);
			  
							/* make sure nzmax increases at least by 1 */
							if (oldnzmax == nLmax) 
								nLmax++;
							mxRealloc( Lvec, nLmax*sizeof(double) );
							mxRealloc( iVec, nLmax*sizeof(mwIndex) );
							mxRealloc( jVec, nLmax*sizeof(mwIndex) );
						}

						x = xe;
						z = ze;
						if ( fabs(z-grz[iz+1])<small ) iz++;
					}
				} else {
					while ( z > zi && z > zr ) {
						iCell = ix*(n_grz-1) + iz;

						ze = grz[iz]>zi ? grz[iz] : zi;
						ze = ze>zr ? ze : zr;
						xe = (ze-b)/m;
						dlx = xe - x;
						dlz = ze - z;
						dl = sqrt( dlx*dlx + dlz*dlz );
			
						/* mexPrintf("x %lf xe %lf      z %lf ze %lf  zi %lf\n", x, xe, z, ze, zi); */

						/* mexPrintf("i = %d, j = %d, dl = %lf, dlx = %lf, dlz = %lf\n", ix, iz, dl, dlx, dlz); */
						iVec[k] = n;
						jVec[k] = iCell;
						Lvec[k] = dl;
						k++;
			
						if (k>=nLmax){
							oldnzmax = nLmax;
							percent_sp += 0.1;
							nLmax = (mwSize)ceil((double)nTx*(double)nSlowness*percent_sp);
			  
							/* make sure nzmax increases at least by 1 */
							if (oldnzmax == nLmax) 
								nLmax++;
							mxRealloc( Lvec, nLmax*sizeof(double) );
							mxRealloc( iVec, nLmax*sizeof(mwIndex) );
							mxRealloc( jVec, nLmax*sizeof(mwIndex) );
						}

						x = xe;
						z = ze;
						if ( fabs(z-grz[iz])<small ) iz--;
					}
				}
		
				ix++;
				x = grx[ix];
			} 
		}
  }
  
  /* ------------------------------------------------------
		 Output variable
     ------------------------------------------------------ */
  nL = k;
  plhs[0] = mxCreateSparse(nTx, nSlowness, nL, mxREAL);
  Lval = mxGetPr( plhs[0] );
  irL  = mxGetIr( plhs[0] );
  jcL  = mxGetJc( plhs[0] );


  k = 0;
  for ( j=0; j<nSlowness; ++j ) {
		jcL[j] = k;

		for ( n=0; n<nL; ++n ) {
			if ( jVec[n] == j ) {
				irL[k] = iVec[n];
				Lval[k] = Lvec[n];
				k++;
			}
		}
  }
  jcL[nSlowness] = k;


  /* ------------------------------------------------------
		 Optional output variables
     ------------------------------------------------------ */
  
  if ( nlhs >= 2 ) {
		plhs[1] = mxCreateDoubleMatrix(n_grx-1, 1, mxREAL);
		gridx = mxGetPr(plhs[1]);
		gridx[0] = grx[0] + 0.5*dx;
		for ( n=1; n<n_grx-1; ++n )
			gridx[n] = gridx[n-1]+dx;
  }
  if ( nlhs >= 3 ) {
		plhs[2] = mxCreateDoubleMatrix(n_grz-1, 1, mxREAL);
		gridz = mxGetPr(plhs[2]);
		gridz[0] = grz[0] + 0.5*dz;
		for ( n=1; n<n_grz-1; ++n )
			gridz[n] = gridz[n-1]+dz;
  }


  mxFree(jVec);
  mxFree(iVec);
  mxFree(Lvec);

  return;
}
