plot for [i=1:50] file=sprintf("output.%d",i) file using 1:10 w l title ""
plot for [i=1:50] file=sprintf("output.%d",i) file using 1:($4+$5) w l title ""
plot for [i=1:50] file=sprintf("scaled.%d",i) file using 1:($12*100) w l title ""
