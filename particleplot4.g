set xlabel 'z(mm)'
set ylabel 'y(mm)'
set zlabel 'x(mm)'
# set zrange [-80:80]
# set yrange [-80:80]
splot 'mirror_simulation4.txt' using 3:2:($4==0 ? $1 : 1/0) pt 7 ps 0.1 title"�Օ��̓��Ŏ��񂾗��q",'mirror_simulation4.txt' using 3:2:(($4==1 && $5==0) ? $1 : 1/0) pt 7 ps 0.1 title"�~���[���Ŏ��񂾗��q",'mirror_simulation4.txt' using 3:2:(($5==1 && $6==0) ? $1 : 1/0) pt 7 ps 0.1 title"�~���[�𔲂��������˂��Ȃ����q",'mirror_simulation4.txt' using 3:2:(($5==1 && $6==1) ? $1 : 1/0) pt 7 ps 0.1 title"���˂������q"
pause -1