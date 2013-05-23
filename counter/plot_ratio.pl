set term postscript color enhanced
set output outputfile

set xlabel "p"
set ylabel "fraction of particles"


plot filename u 1:3 with linespoints title 'Z', filename u 1:4 with linespoints title '{/Symbol s}', filename u 1:5 with linespoints title 'S', filename u 1:2 with linespoints title 'U'
