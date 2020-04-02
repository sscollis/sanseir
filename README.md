# sanseir
Solve a stochastic SEIR model for disease progression

To build with `gfortran` simply:

    make

You will need to modify the `Makefile` if you use a different FORTRAN compiler.
Once you have `sanseir` build, to run simply pass the input file (in FORTRAN 
namelist format) as an argument

    sanseir.exe usa.inp
    
This generates an ensemble of 50 runs with output files called 
`output.*, scaled.*, rates.*`.

These can be plotted using `gnuplot`.  For example, to plot the cummulative deaths

    gnuplot
    plot for [i=1:50] file=sprintf("output.%d",i) file using 1:12 w l lt 'grey' title ""

There are a number of `gnuplot` scripts for making typical plots in the `plot`
directory.

Directory  |  Description
-----------|---------------------------------------------------------------
`Yeo-etal` |  Reproduces results in https://doi.org/10.1101/2020.03.24.20043026 
`test`     |  Additional experimental input files 
`plot`     |  Scripts to make plots using `gnuplot` 


S. Scott Collis
Thu Apr  2 06:42:42 MDT 2020
