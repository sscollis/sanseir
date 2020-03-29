plot for [i=1:50] file=sprintf("output.%d",i) file using 1:10 w l title ""
replot for [i=1:50] file=sprintf("output.%d",i) file using 1:($4+$5) w l title ""
