#
# Plots ensemble of the cumulative number of infected individuals reported
#
set xlabel "Time (days)"
set ylabel "Number of Infected individuals"
set format y "%8.1e"
set mytics 5
set mxtics 5 
set grid
#set yrange [0:1.5e6]
plot for [i=1:1] file=sprintf("output.%d",i) file using 1:($6+0.2*$7+$8+$9) w l lt rgb 'grey10' title "Reported"
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:($6+0.2*$7+$8+$9) w l lt rgb 'grey10' title ""
