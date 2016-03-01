#include <stdio.h>
#include <string.h>

#include "mex.h"

#include <stdint.h>  /* added  */

/**** added for inline void  */
#ifdef _MSC_VER
#define INLINE __forceinline /* use __forceinline (VC++ specific) */
#else
  #define INLINE inline        /* use standard inline */
#endif

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

INLINE void Swap2Bytes(int16_t *x) {   /*inline --> INLINE  */
    *x=(((*x>>8)&0xff) | ((*x&0xff)<<8));
}

INLINE void Swap4Bytes(int32_t *x) {   /*inline --> INLINE  */
    *x=(((*x>>24)&0xff) | ((*x&0xff)<<24) | ((*x>>8)&0xff00) | ((*x&0xff00)<<8));
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

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
	
	  const short NFIELDS_SEG=86;
		short NFIELDS = 0;

		mwSize n, nf, m, ntraces, nfields, number_of_dims, ntraces_one[] = {1, 1};
    const mwSize  *dim_array;
    double *tmp;
    int32_t *traces_no;
    int16_t *fields_no;
    long offset, offset2, ioff;
    
    mxArray **wvalue;
    void    *pdata=NULL;
    
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
    
    FILE *fid;
    char *filename;
    
    if(nrhs<1)
        mexErrMsgTxt("At least one input variables required.");
    else if(nrhs>5)
        mexErrMsgTxt("No more than five input.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    
    int16_t *stmp = (int16_t *)mxMalloc(sizeof(int16_t));
    int32_t *itmp = (int32_t *)mxMalloc(sizeof(int32_t));
    double *dtmp = (double *)mxMalloc(sizeof(double));
    float *ftmp = (float *)mxMalloc(sizeof(float));

    int bt = IsBigEndian();

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
		char msg[100];
    switch( *stmp ) {
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
					sprintf(msg, "Data type (%hd) not compliant to SEG-Y standard", *stmp);
					mexErrMsgTxt(msg);
    }
    
    // get Number of 3200-byte, Extended Textual File Header records
    if ( fseek(fid, 3502, SEEK_SET) == -1 ) {
        fclose(fid);
        mexErrMsgTxt("Problem with the size of the SEG-Y file.");
    }
    fread(stmp, 2, 1, fid);
    if ( bt == false ) {
        Swap2Bytes( stmp );
    }
    int nextended = *stmp;
    


    // arg 2: vector of trace number to read
    //
    ntraces = 0;
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
				
				for ( n = 0; n<NFIELDS; ++n ) {
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
				number_of_dims = mxGetNumberOfDimensions(prhs[4]);
				if ( number_of_dims != 2 ) {
					mexErrMsgTxt("word length must be a rank 2 matrix.");
				}
				dim_array = mxGetDimensions(prhs[4]);
				if( dim_array[0] != 1 && dim_array[1] != 1 ) {
					mexErrMsgTxt("word length must be a vector.");
				}
				nfields = dim_array[0]*dim_array[1];
				if ( nfields > NFIELDS ) {
					mexErrMsgTxt("Number of fields in word length larger than in dictionnary.");
				}
				tmp = (double*)mxGetPr(prhs[4]);
				for (n=0; n<nfields; ++n) {
					word_length[n] = (short) lround( tmp[n] );
				}
				
			} else {
				mexErrMsgTxt("Arg 5 must be of type double.");
			}
			
		} else {

			// assign default dictionnary
			
			NFIELDS = NFIELDS_SEG;
			fnames = (char **)mxMalloc(sizeof(char*)*NFIELDS);
			word_length = (short *)mxMalloc(sizeof(short)*NFIELDS);
			for ( n = 0; n<NFIELDS; ++n ) {
				fnames[n] = (char *)mxMalloc(sizeof(char)*(1+strlen(fnames_seg[n])));
				strcpy(fnames[n], fnames_seg[n]);
				word_length[n] = word_length_seg[n];
			}
		}
		    
    // arg 3: string trace header word number to read
    //
    
    if ( nrhs>=3 ) {
        if (mxIsDouble(prhs[2])) {
            number_of_dims = mxGetNumberOfDimensions(prhs[2]);
            if ( number_of_dims != 2 ){
                mexErrMsgTxt("fields must be a rank 2 matrix.");
            }
            dim_array = mxGetDimensions(prhs[2]);
            if( dim_array[0] != 1 && dim_array[1] != 1 ) {
                mexErrMsgTxt("fields: must be a vector.");
            }
            nfields = dim_array[0]*dim_array[1];
            if ( nfields > NFIELDS ) {
                mexErrMsgTxt("Number of fields larger than allowed by SEG standard.");
            }
            tmp = (double*)mxGetPr(prhs[2]);
            fields_no = (int16_t*)mxMalloc(nfields*sizeof(int16_t));
            for (n=0; n<nfields; ++n)
                fields_no[n] = (int16_t) lround( tmp[n] ) - 1;  // indices start at 0 now
        }
        else if (mxIsChar(prhs[2])) {
            /* copy the string data from prhs[2] into a C string input_ buf.    */
            char *input_buf = mxArrayToString(prhs[2]);
            nfields = 0;
            char *p;
            fields_no = (int16_t*)mxMalloc(NFIELDS*sizeof(int16_t));
            p = strtok(input_buf, " ,");
            while (p != NULL) {
                int found = 0;
                for (n=0; n<NFIELDS; ++n) {
                    if ( strcmp( p, fnames[n] ) == 0 ) {
                        found = 1;
                        fields_no[ nfields ] = n;
                        break;
                    }
                }
                if ( found == 0 ) {
                    char message[80];
                    sprintf(message, "Header field \"%s\" not defined", p);
                    mexErrMsgTxt(message);
                }
                nfields++;
                p = strtok(NULL, " ,");
            }
        } else {
            mexErrMsgTxt("Arg 3 must be of type string or double.");
        }
    } else {
        nfields = NFIELDS;
        fields_no = (int16_t*)mxMalloc(nfields*sizeof(int16_t));
        for (n=0; n<nfields; ++n)
            fields_no[n] = n;
    }

    
    // get the right field names
    //
    char **field_names;
    field_names = mxCalloc(nfields, sizeof(*field_names));
    for (n=0; n<nfields; ++n)
        field_names[n] = mxMalloc( (1+strlen(fnames[fields_no[n]]))*sizeof(char) );
    for (n=0; n<nfields; ++n)
        strcpy( field_names[n], fnames[fields_no[n]] );
    
    /* create a 1 x 1 struct matrix for output  */
    plhs[0] = mxCreateStructMatrix(1, 1, nfields, (const char **)field_names);
    
    ntraces_one[0] = ntraces;
    
    // create the vectors to hold the data
    
    wvalue = (mxArray **)mxMalloc(nfields*sizeof(mxArray *));
    for ( nf=0; nf<nfields; ++nf ) {
        switch ( word_length[ fields_no[nf] ] ) {
            case 2:
                wvalue[nf] = mxCreateNumericArray(2, ntraces_one, mxINT16_CLASS, mxREAL);
                break;
            case 4:
                wvalue[nf] = mxCreateNumericArray(2, ntraces_one, mxINT32_CLASS, mxREAL);
                break;
				    case 5:
							  wvalue[nf] = mxCreateNumericArray(2, ntraces_one, mxSINGLE_CLASS, mxREAL);
                break;
            case 6:
                wvalue[nf] = mxCreateNumericArray(2, ntraces_one, mxDOUBLE_CLASS, mxREAL);
                break;
            default:
                mexErrMsgTxt("This should never happen!");
        }
    }
    
    
    for ( n=0; n<ntraces; ++n ) {
        
        offset = 3600L + nextended*3200 + (traces_no[n]*(240+(bytesPerSample*nsamples)));
        
        for ( nf=0; nf<nfields; ++nf ) {
            
            offset2 = offset;
            for (m=0; m<fields_no[nf]; ++m) {
							  ioff = word_length[m]==5 ? 4 : word_length[m];  // we use 5 for 4-byte ibm float
								offset2 += ioff;
						}
            
            if ( fseek(fid, offset2, SEEK_SET) == -1 ) {
                fclose(fid);
                mexErrMsgTxt("Problem with the size of the SEG-Y file.");
            }
            
            switch ( word_length[ fields_no[nf] ] ) {
                case 2:
                    fread(stmp, 2, 1, fid);
                    
                    if ( bt == false ) {
                        Swap2Bytes( stmp );
                    }
                    
                    pdata = (int16_t *)mxGetData( wvalue[nf] );
                    
                    memcpy(((int16_t*)pdata)+n, stmp, 2);
                    
                    break;
                case 4:
                    fread(itmp, 4, 1, fid);
                    
                    if ( bt == false ) {
                        Swap4Bytes( itmp );
                    }
                    
                    pdata = (int32_t *)mxGetData( wvalue[nf] );
                    
                    memcpy(((int32_t *)pdata)+n, itmp, 4);
                    
                    break;
						    case 5:
									  fread(itmp, 4, 1, fid);
										ibm2float(itmp, (int32_t *)ftmp, 1, bt);

										pdata =(float *) mxGetData( wvalue[nf] );

										memcpy(((float *)pdata)+n, ftmp, sizeof(float));

										break;
                case 6:
                    fread(itmp, 4, 1, fid);
                    fread(stmp, 2, 1, fid);
                    
                    if ( bt == false ) {
                        Swap4Bytes( itmp );
                        Swap2Bytes( stmp );
                    }
                    
                    *dtmp = *itmp * pow( 10.0, *stmp );
                    
                    pdata =(double *) mxGetData( wvalue[nf] );
                    
                    memcpy(((double *)pdata)+n, dtmp, sizeof(double));
                    
                    break;
                default:
                    mexErrMsgTxt("This should never happen!");
            }
        }
    }
    
    //  finally assign the vectors to the fields
    for ( nf=0; nf<nfields; ++nf ) {
        mxSetField(plhs[0], 0, field_names[nf], wvalue[nf]);
    }
    
    fclose(fid);

    return;
}
