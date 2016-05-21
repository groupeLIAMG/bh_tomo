
#include <stdint.h>
#include <string.h>

#include "mex.h"



/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

    FILE *fid;
    char *filename;

    const char *fnames[] = {  // follows CWP/SU naming convention for bytes 3201-3260
    "jobid",  // int
    "lino",   // int
    "reno",   // int
    "ntrpr",  // short
    "nart",   // short
    "hdt",    // unsigned short
    "dto",    // unsigned short
    "hns",    // unsigned short
    "nso",    // unsigned short
    "format", // short
    
    "fold",   // short
    "tsort",  // short
    "vscode", // short
    "hsfs",   // short
    "hsfe",   // short
    "hslen",  // short
    "hstyp",  // short
    "schn",   // short
    "hstas",  // short
    "hstae",  // short
    
    "htatyp", // short
    "hcorr",  // short
    "bgrcv",  // short
    "rcvm",   // short
    "mfeet",  // short
    "polyt",  // short
    "vpol",   // short
    
    // bytes after 3500
    "rev",    // unsigned short
    "fixl",   // short
    "extfh",  // short
    };

    const short word_length[] = {
        4, 4, 4, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2
    };
    
    const short signed_word[] = {
        1, 1, 1, 1, 1, 0, 0, 0, 0, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1
    };
    

    if(nrhs!=2)
      mexErrMsgTxt("Two input variables required.");
    
    else if(nlhs > 0)
      mexErrMsgTxt("Too many output arguments.");

    // arg 1: filename
    //
    if(!mxIsChar(prhs[0]))
      mexErrMsgTxt("Arg 1 must be of type char.");
    filename = mxArrayToString(prhs[0]);
    
    if(!mxIsStruct(prhs[1]))
		mexErrMsgTxt("header must be a structure.");
    
    //  open file
    fid = fopen(filename, "a+");
    if ( fid==NULL ) {
        mexErrMsgTxt("Cannot open SEG-Y file.");
    }

    
    long offset = 3200L;
    if ( fseek(fid, offset, SEEK_SET) == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem writing to SEG-Y file.");
    }
    for ( size_t nf=0; nf<27; ++nf ) {
        
        switch ( word_length[ nf ] ) {
            case 2:
                if ( signed_word[nf] ) {
                    int16_t *val = (int16_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    if ( fwrite(val,2,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file.");
                    }
                } else {
                    uint16_t *val = (uint16_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    if ( fwrite(val,2,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file.");
                    }
                }
                break;
            case 4:
                if ( signed_word[nf] ) {
                    int32_t *val = (int32_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    if ( fwrite(val,4,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file.");
                    }
                } else {
                    uint32_t *val = (uint32_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    if ( fwrite(val,4,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file.");
                    }
                }
                
                break;
            default:
                mexErrMsgTxt("This should never happen!");
        }
    }
    int16_t *zero = (int16_t *)mxMalloc(sizeof(int16_t));
    zero[0] = 0;
    for (size_t n=3260; n<3500; n+=2) {
        if ( fwrite(zero,2,1,fid)!=1) {
            fclose(fid);
            mexErrMsgTxt("Problem writing to SEG-Y file.");
        }
    }
        
    // set SEG Y Format revision number
    
    uint16_t *uval = (uint16_t *)( mxGetData( mxGetField(prhs[1], 0, "rev") ) );
    if ( fwrite(uval,2,1,fid)!=1) {
        fclose(fid);
        mexErrMsgTxt("Problem writing to SEG-Y file.");
    }
    
    // set Fixed length trace flag
    int16_t *val = (int16_t *)( mxGetData( mxGetField(prhs[1], 0, "fixl") ) );
    if ( fwrite(val,2,1,fid)!=1) {
        fclose(fid);
        mexErrMsgTxt("Problem writing to SEG-Y file.");
    }
    
    // set Number of 3200-byte, Extended Textual File Header records
    val = (int16_t *)( mxGetData( mxGetField(prhs[1], 0, "extfh") ) );
    if ( fwrite(val,2,1,fid)!=1) {
        fclose(fid);
        mexErrMsgTxt("Problem writing to SEG-Y file.");
    }
    
    for (size_t n=3506; n<3600; n+=2) {
        if ( fwrite(zero,2,1,fid)!=1) {
            fclose(fid);
            mexErrMsgTxt("Problem writing to SEG-Y file.");
        }
    }
    
    fclose(fid);
    mxFree(zero);
    
    return;
}