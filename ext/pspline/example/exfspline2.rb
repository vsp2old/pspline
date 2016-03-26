#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# complex FFT data points"

XY = [ [0, 1],[0, 1],[0, 1],[1, 0,],[1, 0,],[1, 0],[1, 0],[0, 1],[0, 1],[0, 1] ]

fs = XY.fft_complex_forward
fs.real.each {|x| printf("% .2f ", x) }
puts
fs.image.each {|y| printf("% .2f ", y) }
puts

printf "Real :"
-5.upto(5) {|t|
	w = fs[t]
	printf "% .2f ", w[0]
}
puts
printf "Image:"
-5.upto(5) {|t|
	w = fs[t]
	printf "% .2f ", w[1]
}
puts

bs = fs.spline order:3
vv = fs.axis

s = bs.plot(vv, 10)

puts "# Inverse data"

st = fs.inverse
printf "["
st.each {|x| printf "% .1f ", x[0]}
printf "]\n["
st.each {|x| printf "% .1f ", x[1]}
puts "]"

# Draw Graph
require "gnuplot"
 
Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'sec2(x)'
		plot.ylabel 'Y'
		plot.xlabel 'X'
 		x = vv.map {|v| v[0]}
		y = vv.map {|v| v[1][0]}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.title = "Real"
		end
		z = vv.map {|v| v[1][1]}
		plot.data << Gnuplot::DataSet.new( [x, z] ) do |ds|
			ds.with = "lines"
			ds.title = "Image"
		end
		p = vv.map {|v| a = v[1][0]; b = v[1][1]; Math.sqrt(a*a+b*b)}
		plot.data << Gnuplot::DataSet.new( [x, p] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.title = "Power"
		end
	end
end
