#
# Plots ensemble of the number of infected individuals,
##
set xlabel "Time (days)"
set ylabel "Number of Infected"
set format y "%8.1e"
set mytics 5
set mxtics 5 
set grid
#set yrange [0:1.5e6]
plot for [i=1:1] file=sprintf("output.%d",i) file using 1:4 w l lt rgb 'grey10' title "Ih"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:5 w l lt rgb 'grey50' title "Ic"
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:4 w l lt rgb 'grey10' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:5 w l lt rgb 'grey50' title ""
