#include "mex.h"
#include "string.h"



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
    const short nfields=30;
    mwSize n, nf, m;
    const mwSize oneone[] = {1, 1};
    double *tmp;
    long offset;
    
    mxArray *wvalue;
    void  *pdata=NULL;
    
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
    
    FILE *fid;
    char *filename;
    
    if(nrhs!=1)
      mexErrMsgTxt("One input variables required.");
    
    else if(nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");

    // arg 1: filename
    //
    if(!mxIsChar(prhs[0]))
        mexErrMsgTxt("Arg 1 must be of type char.");
    filename = mxArrayToString(prhs[0]);
    
    
    int16_t *stmp = (int16_t *)mxMalloc(sizeof(int16_t));
    int32_t *itmp = (int32_t *)mxMalloc(sizeof(int32_t));
    double *dtmp = (double *)mxMalloc(sizeof(double));

    int bt = IsBigEndian();

    //  open file
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        mexErrMsgTxt("Cannot open SEG-Y file.");
    }
    

    /* create a 1x1 struct matrix for output  */
    plhs[0] = mxCreateStructMatrix(1, 1, nfields, (const char **)fnames);

    
    offset = 3200L;
    for ( nf=0; nf<27; ++nf ) {
//        if ( word_length[nf] == 240 ||  // unassigned bytes 3261-3500
//                word_length[nf] == 94 ) {  // unassigned bytes 3507-3600
//            offset += word_length[nf];
//            continue;
//        }
        
        if ( fseek(fid, offset, SEEK_SET) == -1 ) {
            fclose(fid);
            mexErrMsgTxt("Problem with the size of the SEG-Y file.");
        }
        
        switch ( word_length[ nf ] ) {
            case 2:
                fread(stmp, 2, 1, fid);
                if ( bt == false ) {
                    Swap2Bytes( stmp );
                }
                if ( signed_word[nf] )
                    wvalue = mxCreateNumericArray(2, oneone, mxINT16_CLASS, mxREAL);
                else
                    wvalue = mxCreateNumericArray(2, oneone, mxUINT16_CLASS, mxREAL);
                
                pdata = mxGetData( wvalue );
                
                memcpy(pdata, stmp, 2);
                
                break;
            case 4:
                fread(itmp, 4, 1, fid);
                if ( bt == false ) {
                    Swap4Bytes( itmp );
                }
                
                if ( signed_word[nf] )
                    wvalue = mxCreateNumericArray(2, oneone, mxINT32_CLASS, mxREAL);
                else
                    wvalue = mxCreateNumericArray(2, oneone, mxUINT32_CLASS, mxREAL);
                
                pdata = mxGetData( wvalue );
                
                memcpy(pdata, itmp, 4);
                
                break;
            default:
                mexErrMsgTxt("This should never happen!");
        }
        
        mxSetField(plhs[0], 0, fnames[ nf ], wvalue);
        
        offset += word_length[nf];
    }
    
    // get SEG Y Format revision number
    
    if ( fseek(fid, 3500L, SEEK_SET) == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
    }

    fread(stmp, 2, 1, fid);
    if ( bt == false ) {
        Swap2Bytes( stmp );
    }
    wvalue = mxCreateNumericArray(2, oneone, mxUINT16_CLASS, mxREAL);
    pdata = mxGetData( wvalue );
    memcpy(pdata, stmp, 2);
    mxSetField(plhs[0], 0, "rev", wvalue);
    
    
    // get Fixed length trace flag
    
    if ( fseek(fid, 3502L, SEEK_SET) == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
    }

    fread(stmp, 2, 1, fid);
    if ( bt == false ) {
        Swap2Bytes( stmp );
    }
    wvalue = mxCreateNumericArray(2, oneone, mxINT16_CLASS, mxREAL);
    pdata = mxGetData( wvalue );
    memcpy(pdata, stmp, 2);
    mxSetField(plhs[0], 0, "fixl", wvalue);
    
    
    // get Number of 3200-byte, Extended Textual File Header records
    
    if ( fseek(fid, 3504L, SEEK_SET) == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
    }

    fread(stmp, 2, 1, fid);
    if ( bt == false ) {
        Swap2Bytes( stmp );
    }
    wvalue = mxCreateNumericArray(2, oneone, mxINT16_CLASS, mxREAL);
    pdata = mxGetData( wvalue );
    memcpy(pdata, stmp, 2);
    mxSetField(plhs[0], 0, "extfh", wvalue);

    
    fclose(fid);
    return;
}
