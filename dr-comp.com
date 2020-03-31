#
# Plots ensemble of death rate with different testing
#
set xlabel "Time (days)"
set ylabel "Death Rate"
set format y "%8.1e"
set mytics 5
set mxtics 5 
set grid
set key top left
plot for [i=1:1] file=sprintf("run1/rates.%d",i) file using 1:12 w l lt rgb 'grey50' lw 2 title "Limited testing"
replot for [i=2:50] file=sprintf("run1/rates.%d",i) file using 1:12 w l lt rgb 'grey50' lw 1 title ""
replot for [i=1:1] file=sprintf("run2/rates.%d",i) file using 1:12 w l lt rgb 'grey10' lw 2 title "Increased testing"
replot for [i=2:50] file=sprintf("run2/rates.%d",i) file using 1:12 w l lt rgb 'grey10' lw 1 title ""
