set terminal pdfcairo enhanced font "Times-New-Roman,25" size 16cm,12cm
set output "../fig/figS2.pdf"

set key font "Times-New-Roman,20"
set key top left


# 軸設定など
set xrange[0:600000]
set xlabel "{/=20  }"
set ylabel "{/=20  }"

set tics font "Times New Roman,25"

plot "f_vwall1.0e-06.txt" w l notitle

unset output