set xlabel 'z(mm)'
set ylabel 'y(mm)'
set zlabel 'x(mm)'
# set zrange [-80:80]
# set yrange [-80:80]
splot 'mirror_simulation4.txt' using 3:2:($4==0 ? $1 : 1/0) pt 7 ps 0.1 title"遮蔽体内で死んだ粒子",'mirror_simulation4.txt' using 3:2:(($4==1 && $5==0) ? $1 : 1/0) pt 7 ps 0.1 title"ミラー内で死んだ粒子",'mirror_simulation4.txt' using 3:2:(($5==1 && $6==0) ? $1 : 1/0) pt 7 ps 0.1 title"ミラーを抜けたが反射しない粒子",'mirror_simulation4.txt' using 3:2:(($5==1 && $6==1) ? $1 : 1/0) pt 7 ps 0.1 title"反射した粒子"
pause -1