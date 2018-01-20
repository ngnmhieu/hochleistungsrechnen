# Plotten von Allem in einer Datei weil Leute das so wollen scheinbar
#====================================
set title 'Weak Scaling Jacobi vs Gauss-Seidel mit Angepasster Problemgroesse'
set yrange [0:1000]
set xrange [0:100]
set logscale x
set pointsize 1
set grid
set xlabel 'Prozessanzahl'
set ylabel 'Berechnungszeit in Sekunden'
plot 'WEAK_SCALING_JA.dat' using 1:4 title "Jacobi" with points, \
'WEAK_SCALING_GS.dat' using 1:4 title "Gauss-Seidel" with points
set term png
set output "weak.png"
replot
set term x11

set title 'Strong Scaling Jacobi vs Gauss-Seidel 960 Interlines 12 Prozesse pro Knoten'
set yrange [0:1000]
set xrange [0:240]
set pointsize 1
set grid
set xlabel 'Anzahl der Prozesse'
set ylabel 'Berechnungszeit in Sekunden'
plot 'STRONG_SCALING_JA.dat' using 1:4 title "Jacobi" with points, \
'STRONG_SCALING_GS.dat' using 1:4 title "Gauss-Seidel" with points
set term png
set output "strong.png"
replot
set term x11

set title 'Communication Jacobi vs Gauss-Seidel 200 Interlines 10 Prozesse'
set yrange [0:1000]
set xrange [0:10]
set pointsize 1
set grid
set xlabel 'Knoten'
set ylabel 'Berechnungszeit in Sekunden'
plot 'COMMUNICATION_A_JA.dat' using 2:4 title "Jacobi" with points, \
'COMMUNICATION_A_GS.dat' using 2:4 title "Gauss-Seidel" with points
set term png
set output "comm.png"
replot
set term x11