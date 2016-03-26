#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

n = 7
j = 5

#
#puts "# Interpolation of the Bessel function"
#
#	Data points
Ad = [ [0.0, 1.0    ],[0.8, 0.84629],[1.6, 0.45540],[2.0, 0.22389],
	   [2.4, 0.00251],[3.2,-0.32019],[4.0,-0.39715], ]
#	Additional data points
C = [ [0.4, 0.96040],[1.2, 0.67113],[2.8,-0.18504],[3.6,-0.39177] ]
Jbn = ARGV[0].to_i
Dp = 10

Bs = Bspline.new(Ad + C, n, j)

vv = []
for i in 0..6
	(p, q) = Ad[i]
	printf "%.2f, % f\n", p, q
	vv.push p
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = Bs.plot(vv, Dp, Jbn) do |u, v|
	printf "% .2f % f\n", u, v[0]
end
# Draw Graph
require "gnuplot"

Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'Bessel'
		plot.ylabel 'Y'
		plot.xlabel 'X'
		x = vv.map {|v| v[0]}
		y = vv.map {|v| v[1]}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end
		y = x.map {|v| Bs.sekibun(v)}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.title = "Integral"
		end
	end
end
