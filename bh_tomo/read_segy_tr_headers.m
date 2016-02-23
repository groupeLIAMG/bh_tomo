% READ_SEGY_TR_HEADERS - read trace header info stored in a SEG-Y file 
%    tr_headers = read_segy_tr_headers(segyfile, traces, fields)
%
% Input:
%    segyfile (mandatory) : name of SEG-Y file
%    traces (optional)  : vector of desired traces within the file, give
%                          empty array to get all traces
%    fields (optional)  : indices or names of trace header word to retreive
%                - indices are given as a vector such as [3 6 8:9]
%                - names are separated by commas such as 'sx,sy,gx,gy'
%
% Output:
%    tr_headers : struct variable containing vectors holding the field
%                 values for the desired traces
%
% Note :  Field names follow the CWP/SU convention for bytes 1-180
%         and are given arbitrarily for bytes 181-240.  For a
%         description of the fields, see Technical Standards
%         (SEG-Y rev 1) at http://www.seg.org 
%
% index   name          bytes in header
% -----   ---------     ---------------
%   1     tracl      -  1-4
%   2     tracr      -  5-8
%   3     fldr       -  9-12
%   4     tracf      -  13-16
%   5     ep         -  17-20
%   6     cdp        -  21-24
%   7     cdpt       -  25-28
%   8     trid       -  29-30
%   9     nvs        -  31-32
%  10     nhs        -  33-34
%
%  11     duse       -  35-36
%  12     offset     -  37-40
%  13     gelev      -  41-44
%  14     selev      -  45-48
%  15     sdepth     -  49-52
%  16     gdel       -  53-56
%  17     sdel       -  57-60
%  18     swdep      -  61-64
%  19     gwdep      -  65-68
%  20     scalel     -  69-70
%
%  21     scalco     -  71-72
%  22     sx         -  73-76
%  23     sy         -  77-80
%  24     gx         -  81-84
%  25     gy         -  85-88
%  26     counit     -  89-90
%  27     wevel      -  91-92
%  28     swevel     -  93-94
%  29     sut        -  95-96
%  30     gut        -  97-98
%
%  31     sstat      -  99-100
%  32     gstat      -  101-102
%  33     tstat      -  103-104
%  34     laga       -  105-106
%  35     lagb       -  107-108
%  36     delrt      -  109-110
%  37     muts       -  111-112
%  38     mute       -  113-114
%  39     ns         -  115-116
%  40     dt         -  117-118
%
%  41     gain       -  119-120
%  42     igc        -  121-122
%  43     igi        -  123-124
%  44     corr       -  125-126
%  45     sfs        -  127-128
%  46     sfe        -  129-130
%  47     slen       -  131-132
%  48     styp       -  133-134
%  49     stas       -  135-136
%  40     stae       -  137-138
%
%  51     tatyp      -  139-140
%  52     afilf      -  141-142
%  53     afils      -  143-144
%  54     nofilf     -  145-146
%  55     nofils     -  147-148
%  56     lcf        -  149-150
%  57     hcf        -  151-152
%  58     lcs        -  153-154
%  59     hcs        -  155-156
%  60     year       -  157-158
%
%  61     day        -  159-160
%  62     hour       -  161-162
%  63     minute     -  163-164
%  64     sec        -  165-166
%  65     timbas     -  167-168
%  66     trwf       -  169-170
%  67     grnors     -  171-172
%  68     grnofr     -  173-174
%  69     grnlof     -  175-176
%  70     gaps       -  177-178
%
%  71     otrav      -  179-180
%  72     xcdp       -  181-184
%  73     ycdp       -  185-188
%  74     ilineno    -  189-192
%  75     clineno    -  193-196
%  76     shotno     -  197-200
%  77     scalsn     -  201-202
%  78     tvmunit    -  203-204
%  79     tdcst      -  205-210
%  80     tdunit     -  211-212
%
%  81     trid       -  213-214
%  82     scalt      -  215-216
%  83     styp       -  217-218
%  84     sdir       -  219-224
%  85     smeas      -  225-230
%  86     smunit     -  231-232
%  87     unass      -  233-240
%
% ---
%
% Bernard Giroux
% INRS-ETE
% 2010-08-04
