#!/usr/bin/env ruby

if ARGV.size < 1
	puts "Usage: ./plot.rb inputfile1 {inputfile2, ...} "
	exit
end

ARGV.each do |inputfile|
	numbers = inputfile.scan( /(\d+(\.\d+)*)/ )
	epsilon = numbers[2][0]
	pressure = numbers[3][0]
	outputname = "e#{epsilon}_p#{pressure}"

	puts "Calculating diffraction data.... for #{inputfile}"
	system "./sofk_2d #{inputfile} #{outputname}.dat"

	puts "Plotting data...."
	system "gnuplot -e \"filename='#{outputname}.dat'; outputfile='#{outputname}.png'\" plotdata.gnuplot"

	puts "Opening output..."
	system "open #{outputname}.png"
end
