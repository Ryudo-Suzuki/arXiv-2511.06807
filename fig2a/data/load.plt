set terminal pdfcairo enhanced font "Times-New-Roman,25" size 16cm,12cm
set output "../fig/fig2a.pdf"

set key font "Times-New-Roman,40"
set key left top offset -5,0


set size 1,1
set origin 0,0
# パラメータ
r = 0.4
ep = 0.01
sqrt3 = sqrt(3)
muc = 2-sqrt3

# 関数定義
z_star(x) = (x**2 + 2 * sqrt3 * x - 1) / (x * (1 + x**2)) / (1 + 1/r * (1 + x**2))
C(x) = (x < muc) ? muc * (x + 2 +sqrt3) /(x * (1+x**2)) / (1 + 1/r * (1 + x**2)) : -2 * (sqrt3 * x**2 - 2 * x - sqrt3 + 1/r * (1 + x**2) * ( x**3 + 3 * sqrt3 * x**2 - 3 * x - sqrt3)) / ((1 + x**2)**2 * (1 + 1/r * (1 + x**2))**2)
z(x,ep) = (z_star(x) + sqrt(z_star(x)**2 + 4 * ep * C(x))) / 2
f(x,ep) = 2 / ep * z(x,ep)

# 軸設定など
set xrange [0.12:0.42]
set yrange [1:10**4]
unset logscale x
set logscale y

set tics font "Times New Roman,25"
set xlabel "{/=30  }"
set ylabel "{/=30  }"

set format x "%g"
set format y "10^{%T}"

# プロット
plot  \
    "line.txt" w vectors nohead dt (3,3) lw 3 lc "gray" notitle,\
    "k1e+2.txt" u 1:($2+3) w p pt 7 ps 1.5 lc 2 title " ",\
    "k1e+3.txt" u 1:($2+3) w p pt 9 ps 1.5 lc 3 title " " ,\
    "k1e+4.txt" u 1:($2+3) w p pt 13 ps 1.5 lc 4 title " " ,\
    f(x,0.01) w l lw 2 lc 2 notitle,\
    f(x,0.001) w l lw 2 lc 3 notitle,\
    f(x,0.0001) w l lw 2 lc 4 notitle,\
    "rigid.txt" u 1:($2+3) w l lw 3 lc 1 notitle

unset output