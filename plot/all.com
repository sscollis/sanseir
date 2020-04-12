#
# Plots ensemble of the number of infected individuals, actuals and reported
#
set xlabel "Time (days)"
set ylabel "Individuals"
set format y "%8.1e"
set mytics 5
set mxtics 5 
set grid
#set yrange [0:1.5e6]
  plot for [i=1:1] file=sprintf("output.%d",i) file using 1:2 w l lt rgb 'grey10' title "S"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:3 w l lt rgb 'grey50' title "E"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:4 w l lt rgb 'grey50' title "Ih"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:5 w l lt rgb 'grey50' title "Ic"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:6 w l lt rgb 'grey50' title "Rh"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:7 w l lt rgb 'grey50' title "Rc"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:8 w l lt rgb 'grey50' title "Dh"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:8 w l lt rgb 'grey50' title "Dc"

replot for [i=2:50] file=sprintf("output.%d",i) file using 1:2 w l lt rgb 'grey10' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:3 w l lt rgb 'grey50' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:4 w l lt rgb 'grey50' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:5 w l lt rgb 'grey50' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:6 w l lt rgb 'grey50' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:7 w l lt rgb 'grey50' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:8 w l lt rgb 'grey50' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:8 w l lt rgb 'grey50' title ""
