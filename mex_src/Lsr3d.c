/*
 * Matrix L from straight rays, in 3D
 *
 * Lsr3d
 *
 *  Created by Bernard Giroux on 13-01-24
 *  Copyright 2013 Bernard Giroux.
 *
 */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <math.h>

#include "mex.h"

#define sign(x) (x > 0) - (x < 0)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double        dx, dy, dz, x1, y1, z1, x2, y2, z2, b_y, m_y, b_z, m_z;
	double        dl, dlx, dly, dlz, x, y, z;
	double        percent_sp, xe, ye, ze, yi, zi;
	const double  small=1.e-10;
	double        *Tx, *Rx, *grx, *gry, *grz, *Lval, *Lvec;
	double        *gridx, *gridy, *gridz;
	double        d, l, m, n;
	mwSize        nTx, nRx, n_grx, n_gry, n_grz, nSlowness;
	mwSize        number_of_dims, nLmax;
	const mwSize  *dim_array;
	mwIndex       *irL, *jcL, *iVec, *jVec;
	mwSize        nt, ix, iy, iz, iCell, i, j, k, nL, oldnzmax;
	int           up_y, up_z;
	double        sx, sy, sz, dtmp;
	
	
	/* Check for proper number of arguments */
	
	if (nrhs != 5) {
		mexErrMsgTxt("LSR3D requires five input arguments.");
	} else if (nlhs > 4) {
		mexErrMsgTxt("LSR3D has a maximum of four output argument.");
	}
	
	/* ------------------------------------------------------
	 * Tx
	 * ------------------------------------------------------ */
	
	if (!(mxIsDouble(prhs[0]))) {
		mexErrMsgTxt("Tx must be double precision.");
	}
	number_of_dims = mxGetNumberOfDimensions(prhs[0]);
	if ( number_of_dims != 2 ){
		mexErrMsgTxt("Tx must be a rank 2 matrix.");
	}
	dim_array = mxGetDimensions(prhs[0]);
	nTx = dim_array[0];
	if ( dim_array[1] != 3 ) {
		mexErrMsgTxt("Tx: matrix nTx by 3.");
	}
	Tx = (double*)mxGetPr(prhs[0]);
	
	/* ------------------------------------------------------
	 * Rx
	 * ------------------------------------------------------ */
	
	if (!(mxIsDouble(prhs[1]))) {
		mexErrMsgTxt("Rx must be double precision.");
	}
	number_of_dims = mxGetNumberOfDimensions(prhs[1]);
	if ( number_of_dims != 2 ){
		mexErrMsgTxt("Rx must be a rank 2 matrix.");
	}
	dim_array = mxGetDimensions(prhs[1]);
	nRx = dim_array[0];
	if ( dim_array[1] != 3 ) {
		mexErrMsgTxt("Rx: matrix nRx by 3.");
	}
	Rx = (double*)mxGetPr(prhs[1]);
	
	if ( nTx != nRx ) {
		mexErrMsgTxt("nTx should be equal to nRx.");
	}
	
	/* ------------------------------------------------------
	 * grx
	 * ------------------------------------------------------ */
	
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
	 * gry
	 * ------------------------------------------------------ */
	
	if (!(mxIsDouble(prhs[3]))) {
		mexErrMsgTxt("gry must be double precision.");
	}
	number_of_dims = mxGetNumberOfDimensions(prhs[3]);
	if ( number_of_dims != 2 ){
		mexErrMsgTxt("gry must be a vector.");
	}
	dim_array = mxGetDimensions(prhs[3]);
	if ( dim_array[0] != 1 && dim_array[1] != 1 ) {
		mexErrMsgTxt("gry: vector 1 by n_gry.");
	}
	n_gry = dim_array[1] > dim_array[0] ? dim_array[1] : dim_array[0];
	gry = (double*)mxGetPr(prhs[3]);
	dy = gry[1]-gry[0];
	
	/* ------------------------------------------------------
	 * grz
	 * ------------------------------------------------------ */
	
	if (!(mxIsDouble(prhs[4]))) {
		mexErrMsgTxt("grz must be double precision.");
	}
	number_of_dims = mxGetNumberOfDimensions(prhs[4]);
	if ( number_of_dims != 2 ){
		mexErrMsgTxt("grz must be a vector.");
	}
	dim_array = mxGetDimensions(prhs[4]);
	if ( dim_array[0] != 1 && dim_array[1] != 1 ) {
		mexErrMsgTxt("grz: vector 1 by n_grz.");
	}
	n_grz = dim_array[1] > dim_array[0] ? dim_array[1] : dim_array[0];
	grz = (double*)mxGetPr(prhs[4]);
	dz = grz[1]-grz[0];
	
	
	nSlowness = (n_grx-1)*(n_gry-1)*(n_grz-1);
	nLmax = nTx * n_grx * n_gry * n_grz/2;
	percent_sp = (nLmax*1.0)/(nTx*nSlowness*1.0);
	
	Lvec = mxMalloc( nLmax*sizeof(double) );
	iVec = mxMalloc( nLmax*sizeof(mwIndex) );
	jVec = mxMalloc( nLmax*sizeof(mwIndex) );
	
	k = 0;
	for ( nt=0; nt<nTx; ++nt ) {
		
		x1 = Tx[nt];
		y1 = Tx[nt+nTx];
		z1 = Tx[nt+2*nTx];
		x2 = Rx[nt];
		y2 = Rx[nt+nRx];
		z2 = Rx[nt+2*nRx];
		
		if ( x1>x2 ) {  // on veut x croissant
			dtmp = x1;
			x1 = x2;
			x2 = dtmp;
			dtmp = y1;
			y1 = y2;
			y2 = dtmp;
			dtmp = z1;
			z1 = z2;
			z2 = dtmp;
		}
		
		d = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
		// cosinus directeurs
		l = (x2-x1)/d;
		m = (y2-y1)/d;
		n = (z2-z1)/d;
		
		sx = sign(l);
		sy = sign(m);
		sz = sign(n);
		
		x = x1;
		y = y1;
		z = z1;
		
		for ( ix=0; ix<n_grx-1; ++ix ) if ( x<grx[ix+1] && x>=grx[ix] ) break;
		for ( iy=0; iy<n_gry-1; ++iy ) if ( y<gry[iy+1] && y>=gry[iy] ) break;
		for ( iz=0; iz<n_grz-1; ++iz ) if ( z<grz[iz+1] && z>=grz[iz] ) break;
		
		if ( fabs(l)<small ) {
			if ( fabs(m)<small ) {
				// X & Y constants
				if ( z1>z2 ) {
					dtmp = z1;
					z1 = z2;
					z2 = dtmp;
					z = z1;
					for ( iz=0; iz<n_grz-1; ++iz ) if ( z<grz[iz+1] && z>=grz[iz] ) break;
				}
				while ( z < z2 ) {
					iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
					
					dlz = ( grz[iz+1]<z2 ? grz[iz+1] : z2 ) - z;
					
					iVec[k] = nt;
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
				
			} else if ( fabs(n)<small ) {
				// X & Z constants
				if ( y1>y2 ) {
					dtmp = y1;
					y1 = y2;
					y2 = dtmp;
					y = y1;
					for ( iy=0; iy<n_gry-1; ++iy ) if ( y<gry[iy+1] && y>=gry[iy] ) break;
				}
				while ( y < y2 ) {
					iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
					
					dly = ( gry[iy+1]<y2 ? gry[iy+1] : y2 ) - y;
					
					iVec[k] = nt;
					jVec[k] = iCell;
					Lvec[k] = dly;
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
					
					iy++;
					y = gry[iy];
				}
				
			} else {
				// seul X constant
				if ( y1>y2 ) {
					dtmp = y1;
					y1 = y2;
					y2 = dtmp;
					dtmp = z1;
					z1 = z2;
					z2 = dtmp;
					y = y1;
					for ( iy=0; iy<n_gry-1; ++iy ) if ( y<gry[iy+1] && y>=gry[iy] ) break;
					z = z1;
					for ( iz=0; iz<n_grz-1; ++iz ) if ( z<grz[iz+1] && z>=grz[iz] ) break;
				}
				
				m_z = (z2-z1)/(y2-y1);
				b_z = z2 - m_z*y2;
				up_z = m_z>0;
				
				while ( y < y2 ) {
					
					zi = m_z*gry[iy+1] + b_z;
					
					if ( up_z ) {
						while ( z < zi && z < z2 ) {
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							ze = grz[iz+1]<zi ? grz[iz+1] : zi;
							ze = ze<z2 ? ze : z2;
							ye = (ze-b_z)/m_z;
							dly = ye - y;
							dlz = ze - z;
							dl = sqrt( dly*dly + dlz*dlz );
							
							iVec[k] = nt;
							jVec[k] = iCell;
							Lvec[k] = dl;
							k++;
							
							if (k>=nLmax) {
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
							
							y = ye;
							z = ze;
							if ( fabs(z-grz[iz+1])<small ) iz++;
						}
					} else { // down
						while ( z > zi && z > z2 ) {
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							ze = grz[iz]>zi ? grz[iz] : zi;
							ze = ze>z2 ? ze : z2;
							ye = (ze-b_z)/m_z;
							dly = ye - y;
							dlz = ze - z;
							dl = sqrt( dly*dly + dlz*dlz );
							
							iVec[k] = nt;
							jVec[k] = iCell;
							Lvec[k] = dl;
							k++;
							
							if (k>=nLmax) {
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
							
							y = ye;
							z = ze;
							if ( fabs(z-grz[iz])<small ) iz--;
						}
					}
					iy++;
					y = gry[iy];
				}
			}
		} else {
			if ( fabs(m)<small && fabs(n)<small ) {
				// Y & Z constants
				while ( x < x2 ) {
					iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
					
					dlx = ( grx[ix+1]<x2 ? grx[ix+1] : x2 ) - x;
					
					iVec[k] = nt;
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
				
			} else if ( fabs(m)<small ) {
				// seul Y constant
				m_z = (z2-z1)/(x2-x1);
				b_z = z2 - m_z*x2;
				up_z = m_z>0;
				
				while ( x < x2 ) {
					
					zi = m_z*grx[ix+1] + b_z;
					
					if ( up_z ) {
						while ( z < zi && z < z2 ) {
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							ze = grz[iz+1]<zi ? grz[iz+1] : zi;
							ze = ze<z2 ? ze : z2;
							xe = (ze-b_z)/m_z;
							dlx = xe - x;
							dlz = ze - z;
							dl = sqrt( dlx*dlx + dlz*dlz );
							
							iVec[k] = nt;
							jVec[k] = iCell;
							Lvec[k] = dl;
							k++;
							
							if (k>=nLmax) {
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
					} else { // down
						while ( z > zi && z > z2 ) {
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							ze = grz[iz]>zi ? grz[iz] : zi;
							ze = ze>z2 ? ze : z2;
							xe = (ze-b_z)/m_z;
							dlx = xe - x;
							dlz = ze - z;
							dl = sqrt( dlx*dlx + dlz*dlz );
							
							iVec[k] = nt;
							jVec[k] = iCell;
							Lvec[k] = dl;
							k++;
							
							if (k>=nLmax) {
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
				
			} else if ( fabs(n)<small ) {
				// seul Z constant
				m_y = (y2-y1)/(x2-x1);
				b_y = y2 - m_y*x2;
				up_y = m_y>0;
				
				while ( x < x2 ) {
					
					yi = m_y*grx[ix+1] + b_y;
					
					if ( up_y ) {
						while ( y < yi && y < y2 ) {
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							ye = gry[iy+1]<yi ? gry[iy+1] : yi;
							ye = ye<y2 ? ye : y2;
							xe = (ye-b_y)/m_y;
							dlx = xe - x;
							dly = ye - y;
							dl = sqrt( dlx*dlx + dly*dly );
							
							iVec[k] = nt;
							jVec[k] = iCell;
							Lvec[k] = dl;
							k++;
							
							if (k>=nLmax) {
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
							y = ye;
							if ( fabs(y-gry[iy+1])<small ) iy++;
						}
					} else { // down
						while ( y > yi && y > y2 ) {
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							ye = gry[iy]>yi ? gry[iy] : yi;
							ye = ye>y2 ? ye : y2;
							xe = (ye-b_y)/m_y;
							dlx = xe - x;
							dly = ye - y;
							dl = sqrt( dlx*dlx + dly*dly );
							
							iVec[k] = nt;
							jVec[k] = iCell;
							Lvec[k] = dl;
							k++;
							
							if (k>=nLmax) {
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
							y = ye;
							if ( fabs(y-gry[iy])<small ) iy--;
						}
					}
					ix++;
					x = grx[ix];
				}
				
			} else {
				
				while ( x < x2 ) {
					
					m_y = (y2-y1)/(x2-x1);
					b_y = y2 - m_y*x2;
					up_y = m_y>0;
					
					m_z = (z2-z1)/(x2-x1);
					b_z = z2 - m_z*x2;
					up_z = m_z>0;
					
					yi = m_y*grx[ix+1] + b_y;
					zi = m_z*grx[ix+1] + b_z;
					
					
					while ( (sy*(yi-y))>0 && (sy*(y2-y))>0 &&
							(sz*(zi-z))>0 && (sz*(z2-z))>0 ) {
						
						
						
						if ( up_y ) {
							ye = gry[iy+1]<yi ? gry[iy+1] : yi;
							ye = ye<y2 ? ye : y2;
						} else {
							ye = gry[iy]>yi ? gry[iy] : yi;
							ye = ye>y2 ? ye : y2;
						}
						if ( up_z ) {
							ze = grz[iz+1]<zi ? grz[iz+1] : zi;
							ze = ze<z2 ? ze : z2;
						} else {
							ze = grz[iz]>zi ? grz[iz] : zi;
							ze = ze>z2 ? ze : z2;
						}
						
						if ( (ze-b_z)/m_z < (ye-b_y)/m_y ) {  // we cross z before y
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							xe = (ze-b_z)/m_z;
							ye = m_y*xe + b_y;
							dlx = xe - x;
							dly = ye - y;
							dlz = ze - z;
							
							dl = sqrt( dlx*dlx + dly*dly + dlz*dlz );
							
							iVec[k] = nt;
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
							y = ye;
							z = ze;
							if ( up_z ) {
								if ( fabs(z-grz[iz+1])<small ) iz++;
							} else {
								if ( fabs(z-grz[iz])<small ) iz--;
							}
							
						} else { // we cross y before z
							iCell = (iz*(n_gry-1)+iy)*(n_grx-1) + ix;
							
							xe = (ye-b_y)/m_y;
							ze = m_z*xe + b_z;
							dlx = xe - x;
							dly = ye - y;
							dlz = ze - z;
							
							dl = sqrt( dlx*dlx + dly*dly + dlz*dlz );
							
							iVec[k] = nt;
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
							y = ye;
							z = ze;
							if ( up_y ) {
								if ( fabs(y-gry[iy+1])<small ) iy++;
							} else {
								if ( fabs(y-gry[iy])<small ) iy--;
							}
						}
					}
					ix++;
					x = grx[ix];
				}
			}
		}
	}
	
	/* ------------------------------------------------------
	 * Output variable
	 * ------------------------------------------------------ */
	nL = k;
	plhs[0] = mxCreateSparse(nTx, nSlowness, nL, mxREAL);
	Lval = mxGetPr( plhs[0] );
	irL  = mxGetIr( plhs[0] );
	jcL  = mxGetJc( plhs[0] );
	
	
	k = 0;
	for ( j=0; j<nSlowness; ++j ) {
		jcL[j] = k;
		
		for ( nt=0; nt<nL; ++nt ) {
			if ( jVec[nt] == j ) {
				irL[k] = iVec[nt];
				Lval[k] = Lvec[nt];
				k++;
			}
		}
	}
	jcL[nSlowness] = k;
	
	
	/* ------------------------------------------------------
	 * Optional output variables
	 * ------------------------------------------------------ */
	
	if ( nlhs >= 2 ) {
		plhs[1] = mxCreateDoubleMatrix(n_grx-1, 1, mxREAL);
		gridx = mxGetPr(plhs[1]);
		gridx[0] = grx[0] + 0.5*dx;
		for ( nt=1; nt<n_grx-1; ++nt )
			gridx[nt] = gridx[nt-1]+dx;
	}
	if ( nlhs >= 3 ) {
		plhs[2] = mxCreateDoubleMatrix(n_gry-1, 1, mxREAL);
		gridy = mxGetPr(plhs[2]);
		gridy[0] = gry[0] + 0.5*dy;
		for ( nt=1; nt<n_gry-1; ++nt )
			gridy[nt] = gridy[nt-1]+dy;
	}
	if ( nlhs >= 4 ) {
		plhs[3] = mxCreateDoubleMatrix(n_grz-1, 1, mxREAL);
		gridz = mxGetPr(plhs[3]);
		gridz[0] = grz[0] + 0.5*dz;
		for ( nt=1; nt<n_grz-1; ++nt )
			gridz[nt] = gridz[nt-1]+dz;
	}
	
	
	mxFree(jVec);
	mxFree(iVec);
	mxFree(Lvec);
	
}
