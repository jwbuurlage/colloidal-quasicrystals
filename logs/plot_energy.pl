set term postscript color enhanced
set output outputfile

set xlabel "500 MCS"
set ylabel "E per particle"


plot filename u ($5/1800) with linespoints title 'E'
