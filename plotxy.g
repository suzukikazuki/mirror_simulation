set xlabel 'y(mm)'
set ylabel 'x(mm)'
set key below
plot 'mirror_simulation4.txt' using 2:(($5==1 && $6==0) ? $1 : 1/0) pt 7 ps 0.1 lt 3 title"���˂���",'mirror_simulation4.txt' using 2:(($5==1 && $6==1) ? $1 : 1/0) pt 7 ps 0.1 lt 7 title"����"
pause -1
