# sanseir
Solve a stochastic SEIR model for disease progression

To build with Gfortran simply:

    make
    
and to run

    sanseir.exe usa.inp
    
This generates an ensemble of 50 runs with output files called `output, scaled, rates`.

These can be plotted using Gnuplot.  For example, to plot the cummulative deaths

    plot for [i=1:50] file=sprintf("output.%d",i) file using 1:12 w l lt 'grey' title ""

