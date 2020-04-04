plot for [i=1:1] file=sprintf("output.%d",i) file using 1:($6+0.19*$7+$8+$9) w l lt rgb 'grey10' title "Reported"
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:($6+0.19*$7+$8+$9) w l lt rgb 'grey10' title ""
replot "new.dat" u 1:2 w p t "Infected (data)"
replot "new.dat" u 1:3 w p t "Deaths (data)"
replot for [i=1:1] file=sprintf("output.%d",i) file using 1:12 w l lt rgb 'grey10' title "Deaths"
replot for [i=2:50] file=sprintf("output.%d",i) file using 1:12 w l lt rgb 'grey10' title ""
