set terminal pdfcairo enhanced font "Times-New-Roman,25" size 16cm,12cm
set output "../fig/figS3b.pdf"

set key font "Times-New-Roman,20"
set key top right


# 軸設定など
set xrange[1e-7:1e-3]
set yrange[:5.4]
set xlabel "{/=20  }"
set ylabel "{/=20  }"
set logscale x
set format x "10^{%T}"

set tics font "Times New Roman,25"

plot "vwall.txt" u 1:2 w lp pt 7 ps 2 notitle,\
    4.886751 w l lw 2 title "rigid"

unset output