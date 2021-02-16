load "settings.gnu"
set terminal png size 640, 480
set output "fits.png"

f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x

fit [log10(700):] f1(x) 'out_f.dat' using (log10($1)):(log10($2)) via a1,b1
fit f2(x) 'out_i.dat' using (log10($1)):(log10($2)) via a2,b2
fit f3(x) 'out_k.dat' using (log10($1)):(log10($2)) via a3,b3

plot 'out_f.dat' w p pt 2 lc 'red' title 'MatMul', 10**(f1(log10(x))) title '' dt 2 lw 1.5 lc 'black', 'out_i.dat' w p pt 4 lc 'cyan' title 'ii-jj-kk', 10**(f2(log10(x))) title '' dt 2 lw 1.5 lc 'blue', 'out_k.dat' w p pt 6 lc 'green' title 'kk-jj-ii', 10**(f3(log10(x))) title '' dt 2 lw 1.5 lc 'orange'

