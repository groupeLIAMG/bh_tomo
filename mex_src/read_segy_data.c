#include <float.h>
#include <stdint.h>
#include <string.h>

#include "mex.h"

#ifdef _WIN32
// MS compilers not C99 compliant
long int lround(const double f) {
    long int lx = ( f < 0.0 ) ? f - 0.5 : f + 0.5;
    return lx;
}
#else
#include <math.h>
#endif




/* returns true if on big endian, else false */
int IsBigEndian() { // grabed at http://unixpapa.com/incnote/byteorder.html
    long one = 1;
    return !(*((char *)(&one)));
}

/*
 * Routines Swap2Bytes, Swap4Bytes, SwapFlaot, ibm2float, int32_t2float,
 * int16_t2float, integer1_to_float taken and adapted from CWP/SU
 *
 * Copyright (c) Colorado School of Mines, 2010,
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following
 * conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * -  Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 * -  Neither the name of the Colorado School of Mines nor the names of
 *    its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * Warranty Disclaimer:
 * THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * Export Restriction Disclaimer:
 * We believe that CWP/SU: Seismic Un*x is a low technology product that does
 * not appear on the Department of Commerce CCL list of restricted exports.
 * Accordingly, we believe that our product meets the qualifications of
 * an ECCN (export control classification number) of EAR99 and we believe
 * it fits the qualifications of NRR (no restrictions required), and
 * is thus not subject to export restrictions of any variety.
 *
 *
 */


/******************************************************************************
Authors: Jens Hartmann,   Institut fur Geophysik, Hamburg, Jun 1993
	 John Stockwell, CWP, Colorado School of Mines, Jan 1994
***************************************************************************/
void Swap2Bytes(int16_t *x) {
    *x=(((*x>>8)&0xff) | ((*x&0xff)<<8));
}

void Swap4Bytes(int32_t *x) {
    *x=(((*x>>24)&0xff) | ((*x&0xff)<<24) | ((*x>>8)&0xff00) | ((*x&0xff00)<<8));
}

void SwapFloat(float *x) {
    int32_t *y = (int32_t *)x;
    *y=(((*y>>24)&0xff) | ((*y&0xff)<<24) | ((*y>>8)&0xff00) | ((*y&0xff00)<<8));
}

void ibm2float(int32_t from[], int32_t to[], size_t n, int endian) {
/***********************************************************************
 * ibm2float - convert between 32 bit IBM and IEEE floating numbers
 ************************************************************************
 * Input:
 * from     input vector
 * to       output vector, can be same as input vector
 * endian   byte order =0 little endian (DEC, PC's)
 * =1 other systems
 *************************************************************************
 * Notes:
 * Up to 3 bits lost on IEEE -> IBM
 *
 * IBM -> IEEE may overflow or underflow, taken care of by
 * substituting large number or zero
 *
 * Only integer shifting and masking are used.
 *************************************************************************
 * Credits: CWP: Brian Sumner,  c.1985
 *************************************************************************/
    register int32_t fconv, fmant, t;
    register size_t i;
    
    for (i = 0;i < n; ++i) {
        
        fconv = from[i];
        
        /* if little endian, i.e. endian=0 do this */
        if (endian == 0) fconv = (fconv << 24) | ((fconv >> 24) & 0xff) |
                ((fconv & 0xff00) << 8) | ((fconv & 0xff0000) >> 8);
        
        if (fconv) {
            fmant = 0x00ffffff & fconv;
            /* The next two lines were added by Toralf Foerster */
            /* to trap non-IBM format data i.e. conv=0 data  */
            if (fmant == 0)
                mexWarnMsgTxt("mantissa is zero data may not be in IBM FLOAT Format !");
            t = (int32_t) ((0x7f000000 & fconv) >> 22) - 130;
            while (!(fmant & 0x00800000)) { --t; fmant <<= 1; }
            if (t > 254) fconv = (0x80000000 & fconv) | 0x7f7fffff;
            else if (t <= 0) fconv = 0;
            else fconv =   (0x80000000 & fconv) | (t << 23)
            | (0x007fffff & fmant);
        }
        to[i] = fconv;
    }
    return;
}

void int32_t2float(int32_t from[], float to[], size_t n, int endian)
/****************************************************************************
Author:	J.W. de Bruijn, May 1995
****************************************************************************/
{
	register size_t i;
    
    if (endian == 0) {
        for (i = 0; i < n; ++i) {
            Swap4Bytes(&from[i]);
            to[i] = (float) from[i];
        }
    } else {
        for (i = 0; i < n; ++i) {
            to[i] = (float) from[i];
        }
    }
    return;
}

void int16_t2float(int16_t from[], float to[], size_t n, int endian)
/****************************************************************************
short_to_float - type conversion for additional SEG-Y formats
*****************************************************************************
Author: Delft: J.W. de Bruijn, May 1995
Modified by: Baltic Sea Reasearch Institute: Toralf Foerster, March 1997
****************************************************************************/
{
	register int i;

	if (endian == 0) {
		for (i = n - 1; i >= 0 ; --i) {
			Swap2Bytes(&from[i]);
			to[i] = (float) from[i];
		}
	} else {
		for (i = n - 1; i >= 0 ; --i)
			to[i] = (float) from[i];
	}
    return;
}

 void integer1_to_float(signed char from[], float to[], int n)
/****************************************************************************
integer1_to_float - type conversion for additional SEG-Y formats
*****************************************************************************
Author: John Stockwell,  2005
****************************************************************************/
{
  	while (n--) {
		to[n] = from[n];
	}
}


/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {

    mwSize n, ns, number_of_dims;
    const mwSize  *dim_array;
    double *tmp;
    int32_t *traces_no;

    FILE *fid;
    char *filename;
    
    if(nrhs<1)
        mexErrMsgTxt("At least one input variables required.");
    else if(nrhs>2)
        mexErrMsgTxt("No more than two input.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    

    // arg 1: filename
    //
    if(!mxIsChar(prhs[0]))
        mexErrMsgTxt("Arg 1 must be of type char.");
    filename = mxArrayToString(prhs[0]);
    
    //  open file
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        mexErrMsgTxt("Cannot open SEG-Y file.");
    }
    
    int16_t *stmp = (int16_t *)mxMalloc(sizeof(int16_t));
    int32_t *itmp = (int32_t *)mxMalloc(sizeof(int32_t));
    
    int bt = IsBigEndian();
    
    // read in some info
    // read in number of samples per data trace
    if ( fseek(fid, 3220, SEEK_SET)  == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
    }
    fread(stmp, 2, 1, fid);
    if ( bt == false ) {
        // we have read a big endian word, we're on a little endian machine
        Swap2Bytes( stmp );
    }
    int nsamples = *stmp;
    
    
    // read in data sample format code
    if ( fseek(fid, 3224, SEEK_SET)  == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
    }
    fread(stmp, 2, 1, fid);
    if ( bt == false ) {
        // we have read a big endian word, we're on a little endian machine
        Swap2Bytes( stmp );
    }
    int bytesPerSample;
    short format = *stmp;
    switch( format ) {
        case 1:
        case 2:
        case 5:
            bytesPerSample = 4;
            break;
        case 3:
            bytesPerSample = 2;
            break;
        case 8:
            bytesPerSample = 1;
            break;
        default:
            mexErrMsgTxt("Data type not compliant to SEG-Y standard");
    }
    
//	uint16_t *utmp = (uint16_t *)mxMalloc(sizeof(uint16_t));
//	if ( fseek(fid, 3500, SEEK_SET) == -1 ) {
//        fclose(fid);
//        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
//    }
//    fread(utmp, 2, 1, fid);
//    if ( bt == false ) {
//        Swap2Bytes( utmp );
//    }
//    int rev_no = *utmp;
	
//    fread(utmp, 2, 1, fid);
//    if ( bt == false ) {
//        Swap2Bytes( utmp );
//    }
//    int fixed_length = *utmp;
	
    // get Number of 3200-byte, Extended Textual File Header records
    if ( fseek(fid, 3504, SEEK_SET) == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
    }
    fread(stmp, 2, 1, fid);
    if ( bt == false ) {
        Swap2Bytes( stmp );
    }
    int nextended = *stmp;
    
//	printf("%d  %d  %d\n", rev_no, fixed_length, nextended);
    
    // arg 2: vector of trace number to read
    //
    int ntraces = 0;
    if ( nrhs>=2 ) {
        if (!mxIsEmpty(prhs[1])) {
            if (!mxIsDouble(prhs[1]))
                mexErrMsgTxt("Arg 2 must be of type double.");
            number_of_dims = mxGetNumberOfDimensions(prhs[1]);
            if ( number_of_dims != 2 ){
                mexErrMsgTxt("Rx must be a rank 2 matrix.");
            }
            dim_array = mxGetDimensions(prhs[1]);
            if ( dim_array[0] != 1 && dim_array[1] != 1 ) {
                mexErrMsgTxt("traces: must be a vector.");
            }
            ntraces = dim_array[0]*dim_array[1];
            tmp = (double*)mxGetPr(prhs[1]);
            traces_no = (int32_t*)mxMalloc(ntraces*sizeof(int32_t));
            for (n=0; n<ntraces; ++n)
                traces_no[n] = (int32_t) lround( tmp[n] ) - 1;  // indices start at 0 now
        }
    } 
    if ( ntraces == 0 ) {
        // we must read for all traces
        
        fseek(fid, 0L, SEEK_END);
        long filesize = ftell(fid);
        fseek(fid, 0L, SEEK_SET);
    
        filesize -= 3600L + nextended*3200L;
        ntraces = filesize / ( 240 + bytesPerSample*nsamples );
        traces_no = (int32_t*)mxMalloc(ntraces*sizeof(int32_t));
        for (n=0; n<ntraces; ++n)
            traces_no[n] = n;
    }

    plhs[0] = mxCreateNumericMatrix(nsamples, ntraces, mxSINGLE_CLASS, mxREAL);
    float *pdata = mxGetData( plhs[0] );
    
    void *ptr;
    switch ( format ) {
        case 1:
        case 2:
            ptr = (int32_t *) mxMalloc(nsamples*sizeof(int32_t));
            break;
        case 3:
            ptr = (int16_t *) mxMalloc(nsamples*sizeof(int16_t));
            break;
        case 5:
            ptr = (float *) mxMalloc(nsamples*sizeof(float));
            break;
        case 8:
            ptr = (char *) mxMalloc(nsamples*sizeof(char));
            break;
        default:
            mexErrMsgTxt("Data format not defined in SEG-Y file.");
    }
    
    float *fltptr;
    if ( format == 5 )
        fltptr = ptr;
    else
        fltptr = mxMalloc(nsamples*sizeof(float));
    
    long offset;
    for ( n=0; n<ntraces; ++n ) {
        
        offset = 3600L + nextended*3200 +
                (traces_no[n]*(240+(bytesPerSample*nsamples))) + 240L;

        if ( fseek(fid, offset, SEEK_SET) == -1 ) {
            fclose(fid);
            mexErrMsgTxt("Problem with the size of the SEG-Y file.");
        }

        fread(ptr, bytesPerSample, nsamples, fid);
        
        switch ( format ) {
            case 1:
                ibm2float(ptr, (int32_t *)fltptr, nsamples, bt);
                break;
            case 2:
                int32_t2float((int32_t *)ptr, fltptr, nsamples, bt);
                break;
            case 3:
                int16_t2float((int16_t *)ptr, fltptr, nsamples, bt);
                break;
            case 5:
                if ( bt == false )
                    for ( ns=0; ns<nsamples; ++ns )
                        SwapFloat(&fltptr[ns]);
                break;
            case 8:
                integer1_to_float((signed char *)ptr, fltptr, nsamples);
                break;
        }
        memcpy(pdata, fltptr, bytesPerSample*nsamples);
        pdata += nsamples;
    }
    
    return;
}
