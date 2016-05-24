
#include <stdint.h>
#include <string.h>

#include "mex.h"

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
    
    int bt = IsBigEndian();
    
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
    
    int16_t i16[1];
    int32_t i32[1];
    
    for ( size_t nf=0; nf<27; ++nf ) {
        
        switch ( word_length[ nf ] ) {
            case 2:
                if ( signed_word[nf] ) {
                    int16_t *val = (int16_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    i16[0] = *val;
                    if ( bt == false ) {
                        Swap2Bytes( i16 );
                    }
                    
                    if ( fwrite(i16,2,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file.");
                    }
                } else {
                    uint16_t *val = (uint16_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    i16[0] = *val;
                    if ( bt == false ) {
                        Swap2Bytes( i16 );
                    }

                    if ( fwrite(i16,2,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file.");
                    }
                }
                break;
            case 4:
                if ( signed_word[nf] ) {
                    int32_t *val = (int32_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    i32[0] = *val;
                    if ( bt == false ) {
                        Swap4Bytes( i32 );
                    }
                    
                    if ( fwrite(i32,4,1,fid)!=1) {
                        fclose(fid);
                        mexErrMsgTxt("Problem writing to SEG-Y file.");
                    }
                } else {
                    uint32_t *val = (uint32_t *)( mxGetData( mxGetField(prhs[1], 0, fnames[nf]) ) );
                    i32[0] = *val;
                    if ( bt == false ) {
                        Swap4Bytes( i32 );
                    }

                    if ( fwrite(i32,4,1,fid)!=1) {
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
    i16[0] = *uval;
    if ( bt == false ) {
        Swap2Bytes( i16 );
    }
    if ( fwrite(i16,2,1,fid)!=1) {
        fclose(fid);
        mexErrMsgTxt("Problem writing to SEG-Y file.");
    }
    
    // set Fixed length trace flag
    int16_t *val = (int16_t *)( mxGetData( mxGetField(prhs[1], 0, "fixl") ) );
    i16[0] = *val;
    if ( bt == false ) {
        Swap2Bytes( i16 );
    }
    if ( fwrite(i16,2,1,fid)!=1) {
        fclose(fid);
        mexErrMsgTxt("Problem writing to SEG-Y file.");
    }
    
    // set Number of 3200-byte, Extended Textual File Header records
    val = (int16_t *)( mxGetData( mxGetField(prhs[1], 0, "extfh") ) );
    i16[0] = *val;
    if ( bt == false ) {
        Swap2Bytes( i16 );
    }
    if ( fwrite(i16,2,1,fid)!=1) {
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