#include <stdint.h>
#include <string.h>

#include "mex.h"

#ifdef _WIN32
// MS compilers not C99 compliant
long int lround(const double x) {
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

void float2ibm(const int32_t from[], int32_t to[], size_t n, int endian)
/**********************************************************************
 float_to_ibm - convert between 32 bit IBM and IEEE floating numbers
***********************************************************************
Input:
from       input vector
n          number of floats in vectors
endian     =0 for little endian machine, =1 for big endian machines
Output:
to         output vector, can be same as input vector
***********************************************************************
Notes:
Up to 3 bits lost on IEEE -> IBM
IBM -> IEEE may overflow or underflow, taken care of by
substituting large number or zero
Only integer shifting and masking are used.
***********************************************************************
Credits:     CWP: Brian Sumner
***********************************************************************/
{
    register int32_t fconv, fmant, t;
    register size_t i;

    for (i=0;i<n;++i) {
        fconv = from[i];
        if (fconv) {
            fmant = (0x007fffff & fconv) | 0x00800000;
            t = (int) ((0x7f800000 & fconv) >> 23) - 126;
            while (t & 0x3) { ++t; fmant >>= 1; }
            fconv = (0x80000000 & fconv) | (((t>>2) + 64) << 24) | fmant;
        }
        if(endian==0)
                fconv = (fconv<<24) | ((fconv>>24)&0xff) |
                        ((fconv&0xff00)<<8) | ((fconv&0xff0000)>>8);

        to[i] = fconv;
    }
    return;
}

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{

    const short NFIELDS_SEG=86;
    short NFIELDS = 0;
    int nextended = 0;  // this is set in write_segy.m
    int bytesPerSample = 4;   // this is set in write_segy.m
    FILE *fid;
    char *filename;

    const char *fnames_seg[] = {  // follows CWP/SU naming convention for bytes 1-180
        "tracl",  // 1
        "tracr",
        "fldr",
        "tracf",
        "ep",
        "cdp",
        "cdpt",
        "trid",
        "nvs",
        "nhs",

        "duse",  // 11
        "offset",
        "gelev",
        "selev",
        "sdepth",
        "gdel",
        "sdel",
        "swdep",
        "gwdep",
        "scalel",

        "scalco",  // 21
        "sx",
        "sy",
        "gx",
        "gy",
        "counit",
        "wevel",
        "swevel",
        "sut",
        "gut",

        "sstat", // 31
        "gstat",
        "tstat",
        "laga",
        "lagb",
        "delrt",
        "muts",
        "mute",
        "ns",
        "dt",

        "gain",  // 41
        "igc",
        "igi",
        "corr",
        "sfs",
        "sfe",
        "slen",
        "styp",
        "stas",
        "stae",

        "tatyp",  // 51
        "afilf",
        "afils",
        "nofilf",
        "nofils",
        "lcf",
        "hcf",
        "lcs",
        "hcs",
        "year",

        "day",  // 61
        "hour",
        "minute",
        "sec",
        "timbas",
        "trwf",
        "grnors",
        "grnofr",
        "grnlof",
        "gaps",

        "otrav",  // 71

        //  names below arbitrarily given
        "xcdp",   // 72 - X coord of ensemble (CDP) position of this trace
        "ycdp",   // 73 - Y coord of ensemble (CDP) position of this trace
        "ilineno",
        "clineno",
        "shotno",   // 76 - shotpoint number
        "scalsn",
        "tvmunit",   // 78 - trace value measurement units
        "tdcst",
        "tdunit",

        "trid",   // 81 - device/trace identifier
        "scalt",
        "styp",   // 83 - source type/orientation
        "sdir",
        "smeas",
        "smunit"   // 86 - source measurement units
    };
    char **fnames;

    const short word_length_seg[] = {
        4, 4, 4, 4, 4, 4, 4, 2, 2, 2,
        2, 4, 4, 4, 4, 4, 4, 4, 4, 2,
        2, 4, 4, 4, 4, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 4, 4, 4, 4, 4, 2, 2, 6, 2,
        2, 2, 2, 6, 6, 2
    };
    short *word_length;

    if(nrhs<3)
      mexErrMsgTxt("At least three input variables required.");

    else if(nlhs > 0)
      mexErrMsgTxt("Too many output arguments.");

    int bt = IsBigEndian();

    // arg 1: filename
    //
    if(!mxIsChar(prhs[0]))
      mexErrMsgTxt("Arg 1 must be of type char.");
    filename = mxArrayToString(prhs[0]);


    // arg 2: trace header data
    //
    if(!mxIsStruct(prhs[1]))
		mexErrMsgTxt("header must be a structure.");

    // arg 3: traces
    //
    float *tr = (float *)mxGetData( prhs[2] );
    size_t nsamples = mxGetM(prhs[2]);
    size_t ntraces = mxGetN(prhs[2]);
    float *traces = mxMalloc(ntraces*nsamples*sizeof(float));
    memcpy(traces, tr, ntraces*nsamples*sizeof(float));
    if ( bt == false ) {
        for ( size_t ns=0; ns<ntraces*nsamples; ++ns )
            SwapFloat(&traces[ns]);
    }

    // arg 4 & 5
    // override standard SEG trace header structure
    if ( nrhs==5 ) {
        // arg 4: dictionnary
        // arg 5: word length (code 5 is for IMB 4-byte float)

        if ( mxIsCell(prhs[3]) ) {
            NFIELDS = mxGetNumberOfElements(prhs[3]);

            fnames = (char **)mxMalloc(sizeof(char*)*NFIELDS);
            word_length = (short *)mxMalloc(sizeof(short)*NFIELDS);

            const mxArray *cell_element_ptr;
            char tmpstr[100];

            for ( size_t n = 0; n<NFIELDS; ++n ) {
                cell_element_ptr = mxGetCell(prhs[3], n);
                if (cell_element_ptr == NULL) {
                    mexErrMsgTxt("\tEmpty Cell\n");
                }
                if (mxIsChar(cell_element_ptr)) {
                    mxGetString(cell_element_ptr, tmpstr, 99);
                    fnames[n] = (char *)mxMalloc(sizeof(char)*(1+strlen(tmpstr)));
                    strcpy(fnames[n], tmpstr);
                } else {
                    mexErrMsgTxt("cell elements must be char.");
                }
            }
        } else {
            mexErrMsgTxt("dictionnary must be of type cell.");
        }

        if (mxIsDouble(prhs[4])) {
            mwSize number_of_dims = mxGetNumberOfDimensions(prhs[4]);
            if ( number_of_dims != 2 ) {
                mexErrMsgTxt("word length must be a rank 2 matrix.");
            }
            const mwSize *dim_array = mxGetDimensions(prhs[4]);
            if( dim_array[0] != 1 && dim_array[1] != 1 ) {
                mexErrMsgTxt("word length must be a vector.");
            }
            mwSize nfields = dim_array[0]*dim_array[1];
            if ( nfields > NFIELDS ) {
                mexErrMsgTxt("Number of fields in word length larger than in dictionnary.");
            }
            double *tmp = (double*)mxGetPr(prhs[4]);
            for (size_t n=0; n<nfields; ++n) {
                word_length[n] = (short) lround( tmp[n] );
            }

        } else {
            mexErrMsgTxt("Arg 4 must be of type double.");
        }

    } else {

        // assign default dictionnary

        NFIELDS = NFIELDS_SEG;
        fnames = (char **)mxMalloc(sizeof(char*)*NFIELDS);
        word_length = (short *)mxMalloc(sizeof(short)*NFIELDS);
        for ( size_t n = 0; n<NFIELDS; ++n ) {
            fnames[n] = (char *)mxMalloc(sizeof(char)*(1+strlen(fnames_seg[n])));
            strcpy(fnames[n], fnames_seg[n]);
            word_length[n] = word_length_seg[n];
        }
    }



    int header_length=0;
    for ( size_t nf=0; nf<NFIELDS; ++nf ) {
        if (word_length[nf] == 5 )
            header_length += 4;
        else
            header_length += word_length[nf];
    }
    if ( header_length>240 ) {
        mexErrMsgTxt("Trace header should be <= 240 bytes.");
    }

    //  open file
    fid = fopen(filename, "a+");
    if ( fid==NULL ) {
        mexErrMsgTxt("Cannot open SEG-Y file.");
    }

    long offset = 3600L + nextended*3200;

    if ( fseek(fid, offset, SEEK_SET) == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem writing to SEG-Y file 0.");
    }

    int8_t *zero = (int8_t *)mxMalloc(sizeof(int8_t));
    zero[0] = 0;

    int16_t **i16v = (int16_t **)mxMalloc(NFIELDS*sizeof(int16_t*));
    int32_t **i32v = (int32_t **)mxMalloc(NFIELDS*sizeof(int32_t*));
    float   **fv = (float **)mxMalloc(NFIELDS*sizeof(float*));
    double  **dv = (double **)mxMalloc(NFIELDS*sizeof(double*));
    for ( size_t nf=0; nf<NFIELDS; ++nf ) {
        
        switch ( word_length[ nf ] ) {
            case 2:
                i16v[nf] = (int16_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                break;
            case 4:
                i32v[nf] = (int32_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                break;
            case 5:
                fv[nf] = (float *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                break;
            case 6:
                dv[nf] = (double *)( mxGetPr( mxGetField(prhs[1], 0, fnames[nf]) ) );
                break;
            default:
                mexErrMsgTxt("This should never happen!");
        }
    }
    
    int16_t i16[1];
    int32_t i32[1];
    double e, m;
    for ( size_t n=0; n<ntraces; ++n ) {
        //mexPrintf("%zd - %ld  --  ", n, ftell(fid));
        // write header
        for ( size_t nf=0; nf<NFIELDS; ++nf ) {
            
            switch ( word_length[ nf ] ) {
                case 2:
                    {
                    i16[0] = i16v[nf][n];
                    if ( bt == false ) {
                        Swap2Bytes( i16 );
                    }
                    if ( fwrite(i16,2,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file 1.");
                    }
                    break;
                    }
                case 4:
                    {
                    i32[0] = i32v[nf][n];
                    //if ( nf==0 ) mexPrintf("%d\n", i32v[nf][n]);
                    if ( bt == false ) {
                        Swap4Bytes( i32 );
                    }
                    if ( fwrite(i32,4,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file 2.");
                    }
                    break;
                    }
                case 5:
                    float2ibm((int32_t *)&(fv[nf][n]), i32, 1, bt);
                    if ( fwrite(i32,4,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file 3.");
                    }
                    break;
                case 6:
                    e = floor(log10(dv[nf][n]))-8;
                    m = dv[nf][n] / pow( 10.0, e);
                    int32_t im[] = { (int32_t)m };
                    int16_t ie[] = { (int16_t)e };
                    if ( bt == false ) {
                        Swap4Bytes( im );
                        Swap2Bytes( ie );
                    }
                    if ( fwrite(im,4,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file 4.");
                    }
                    if ( fwrite(ie,2,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file 5.");
                    }
                    break;
                default:
                    mexErrMsgTxt("This should never happen!");
            }
        }
        // fill remaining of header with 1-byte 0
        for (size_t nn=header_length; nn<240; ++nn) {
            if ( fwrite(zero,1,1,fid)!=1) {
                fclose(fid);
                mexErrMsgTxt("Problem writing to SEG-Y file 6.");
            }
        }

        // write traces
        if ( fwrite(&(traces[n*nsamples]),4,nsamples,fid)!=nsamples) {
            fclose(fid);
            mexErrMsgTxt("Problem writing to SEG-Y file 7.");
        }

    }

    fclose(fid);
    mxFree(zero);
    mxFree(word_length);
    for ( size_t n = 0; n<NFIELDS; ++n ) {
        mxFree(fnames[n]);
    }
    mxFree(fnames);
    mxFree(i16v);
    mxFree(i32v);
    mxFree(fv);
    mxFree(dv);
    mxFree(traces);
    return;
}
