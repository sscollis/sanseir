#
# Plots ensemble of death rate
#
set xlabel "Time (days)"
set ylabel "Death Rate"
set format y "%8.1e"
set mytics 5
set mxtics 5 
set grid
set key top left
plot for [i=1:1] file=sprintf("rates.%d",i) file using 1:12 w l lt rgb 'grey50' title ""
replot for [i=2:50] file=sprintf("rates.%d",i) file using 1:12 w l lt rgb 'grey50' title ""
