# bh_tomo

bh_tomo is an open source borehole georadar data processing and
ray-based 2D and 3D tomography software package.

bh_tomo runs on Matlab 7.0.4 (Release 14 Service Pack 2) or later.

Watch for manual_bh_tomo.pdf and

Giroux, B., Gloaguen E. and Chouteau M.; bh_tomo - A Matlab borehole
georadar 2D tomography package. Computers and Geosciences, 33,
126--137, 2007. doi:10.1016/j.cageo.2006.05.014

for a detailed description.

## Installation

- Copy the files in this directory somewhere in your MATLABPATH or add
the bh_tomo directory to your MATLABPATH environment variable.

- In order to use the S-transform in bh_tomo_amp, you will need to download

www.cora.nwra.com/~stockwel/rgspages/S-Transform/m_files/st.m

save it in your MATLABPATH and **rename it st_bg.m**.  You can turn
verbose off in this file to avoid numerous message printed in your
Matlab command window.

I have included a modified version of the suptitle.m file that allows
to pass optional arguments to the text object.


### MEX files

You might need to compile mex files if not already available for
your platform.  To do so, follow the steps:
```
cd /path/to/bh_tomo/mex_src    ( where the source code is )
mex -O -largeArrayDims Lsr2d.c
mex -O -largeArrayDims Lsr3d.c
mex -O -largeArrayDims Lsr2da.c
mex -O -largeArrayDims ttcr2d.cpp
mex -O -largeArrayDims ttcr2da.cpp
mex -O -largeArrayDims ttcr2daa.cpp
mex -O -largeArrayDims grid3d_mex.cpp
mex -O -largeArrayDims read_segy_b_header.c
mex -O -largeArrayDims read_segy_tr_headers.c
mex -O -largeArrayDims read_segy_data.c
mv *.mex* ../bh_tomo   % (to put the mex files in the main directory of bh_tomo)
```

Invocation of the command mex can be performed in matlab, at the
command prompt.

**Important**: you need a compiler that supports the C++11 standard to
    compile most of the mex functions.  This has to be enabled in your
    mex setup.

## License

bh_tomo is released under the terms of the GNU GENERAL PUBLIC LICENSE.
See LICENSE.txt for details.


Enjoy!

Bernard Giroux
