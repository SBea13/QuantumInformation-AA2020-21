load "settings.gnu"
set terminal png size 640, 480
set output "data.png"

plot 'out_f.dat' w p pt 2 lc 'red' title 'MatMul', 'out_i.dat' w p pt 4 lc 'cyan' title 'ii-jj-kk', 'out_k.dat' w p pt 6 lc 'green' title 'kk-jj-ii'

