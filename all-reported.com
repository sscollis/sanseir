#
# Plots ensemble of the number of infected individuals, actuals and reported
#
set xlabel "Time (days)"
set ylabel "Number of Infected"
set format y "%8.1e"
set mytics 5
set mxtics 5 
set grid
plot for [i=1:1] file=sprintf("output.%d",i) file using 1:10 w l lt rgb 'grey10' title "Reported"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:($4+$5) w l lt rgb 'grey50' title "Actual"
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:10 w l lt rgb 'grey10' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:($4+$5) w l lt rgb 'grey50' title ""
