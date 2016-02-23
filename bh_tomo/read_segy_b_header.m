% READ_SEGY_B_HEADER - read binary header info stored in a SEG-Y file 
%    b_headers = read_segy_b_header(segyfile)
%
% Input:
%    segyfile  : name of SEG-Y file
%
% Output:
%    b_headers : struct variable containing the header field values
%
%
% Note :  Field names follow the CWP/SU convention for bytes 3201-3260
%         and are given arbitrarily for remaining bytes.  For a
%         description of the fields, see Technical Standards
%         (SEG-Y rev 1) at http://www.seg.org 
%
% index   name          bytes in file
% -----   ---------     ---------------
%    1    jobid            3201-3204
%    2    lino             3205-3208
%    3    reno             3209-3212
%    4    ntrpr            3213-3214
%    5    nart             3215-3216
%    6    hdt              3217-3218
%    7    dto              3219-3220
%    8    hns              3221-3222
%    9    nso              3223-3224
%   10    format           3225-3226
%
%   11    fold             3227-3228
%   12    tsort            3229-3230
%   13    vscode           3231-3232
%   14    hsfs             3233-3234
%   15    hsfe             3235-3236
%   16    hslen            3237-3238
%   17    hstyp            3239-3240
%   18    schn             3241-3242
%   19    hstas            3243-3244
%   20    hstae            3245-3246
%
%   21    htatyp           3247-3248
%   22    hcorr            3249-3250
%   23    bgrcv            3251-3252
%   24    rcvm             3253-3254
%   25    mfeet            3255-3256
%   26    polyt            3257-3258
%   27    vpol             3259-3260
%
%   28    rev              3501-3502
%   29    fixl             3503-3504
%   30    extfh            3505-3506
%
% ---
%
% Bernard Giroux
% INRS-ETE
% 2010-08-04
