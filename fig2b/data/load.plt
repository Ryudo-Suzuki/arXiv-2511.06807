set terminal pdfcairo enhanced font "Times-New-Roman,25" size 16cm,12cm
set output "../fig/fig2b.pdf"

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
z_star(k) = (k**2 + 2 * sqrt3 * k - 1) / (k * (1 + k**2)) / (1 + 1/r * (1 + k**2))
C(k) = (k < muc) ? muc * (k + 2 +sqrt3) /(k * (1+k**2)) / (1 + 1/r * (1 + k**2)) : -2 * (sqrt3 * k**2 - 2 * k - sqrt3 + 1/r * (1 + k**2) * ( k**3 + 3 * sqrt3 * k**2 - 3 * k - sqrt3)) / ((1 + k**2)**2 * (1 + 1/r * (1 + k**2))**2)
z(k,x) = (z_star(k) + sqrt(z_star(k)**2 + 4 * x * C(k))) / 2

# 軸設定など
set logscale xy

set tics font "Times New Roman,25"
set xlabel "{/=30  }"
set ylabel "{/=30  }"

set format x "10^{%T}"
set format y "10^{%T}"

plot \
    "mu0.25.txt" u 1:($2+3) w p pt 7 ps 1.5 lc 5 title " ",\
    "muc.txt" u 1:($2+3) w p pt 9 ps 1.5 lc 6 title " " ,\
    "mu0.33.txt" u 1:($2+3) w p pt 13 ps 1.5 lc 7 title " " ,\
    2* x * z(0.25, 1.0/x) lc 5 lw 2 notitle,\
    2* x * z(muc, 1.0/x) lc 6 lw 2 notitle,\
    2* x * z(0.33, 1.0/x) lc 7 lw 2 notitle


unset output