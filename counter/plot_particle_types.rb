#!/usr/bin/env ruby

if ARGV.size < 1
	puts "Usage: ./plot.rb inputfile1 {inputfile2, ...} "
	exit
end

ARGV.each do |inputfile|
	outputname = inputfile.chomp(".txt")

	puts "Plotting data...."
	system "gnuplot -e \"filename='#{outputname}.txt'; outputfile='#{outputname}.eps'\" plot_ratio.pl"

    system "ps2pdf #{outputname}.eps"

	puts "Opening output..."
	system "open #{outputname}.pdf"
end
