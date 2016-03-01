# bh_tomo

bh_tomo is an open source borehole georadar/seismic data processing and
ray-based 2D and 3D tomography software package.

## Installation

Copy the files in the bh_tomo directory somewhere in your MATLABPATH or add
the bh_tomo directory to your MATLABPATH.

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

## Documentation

Watch for manual_bh_tomo.pdf to learn how to use it and look at the following
papers to understand the underlying theory:

- Giroux, B., Gloaguen E. and Chouteau M., 2007.  bh_tomo - A Matlab borehole
georadar 2D tomography package. Computers and Geosciences, 33,
126--137. doi:10.1016/j.cageo.2006.05.014

- Giroux, B., Bouchedda, A. et Chouteau, M., 2009. Assisted traveltime picking
of crosshole GPR data, Geophysics, 74 (4), J35-J48. doi : 10.1190/1.3141002

- Giroux, B. et Chouteau, M., 2010. Quantitative analysis of water content
estimation errors using Ground Penetrating Radar data and a low-loss
approximation, Geophysics 75, WA241-WA249. doi:10.1190/1.3464329

- Giroux B et Gloaguen E, 2012. Geostatistical traveltime tomography in
 elliptically anisotropic media. Geophysical Prospecting, 60, 1133-1149.
 doi : 10.1111/j.1365-2478.2011.01047.x

- Giroux B et Bouchedda A, 2015. Ray-based time-lapse traveltime tomography,
 SEG Technical Program Expanded Abstracts 2015 : pp. 5466-5471.
 doi : 10.1190/segam2015-5815316.1


## License

bh_tomo is released under the terms of the GNU GENERAL PUBLIC LICENSE.
See LICENSE.txt for details.


Enjoy!

Bernard Giroux
