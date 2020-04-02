#
# Plots ensemble of the number of deaths and recovered
#
set xlabel "Time (days)"
set ylabel "Number of individuals"
set format y "%8.1e"
set mytics 5
set mxtics 5 
set grid
set key top left
set auto
set log y
plot for [i=1:1] file=sprintf("output.%d",i) file using 1:12 w l lt rgb 'grey10' title "Deaths"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:11 w l lt rgb 'grey50' title "Recovered"
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:12 w l lt rgb 'grey10' title ""
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:11 w l lt rgb 'grey50' title ""
