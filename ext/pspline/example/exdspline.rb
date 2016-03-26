#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# Interpolation of the sin(x) function"

XY = [ [0.0],[1.0],[2.0],[3.0],[4.0],[5.0],[6.283185] ]

XY.each do |a|
	a.push Math.sin(a[0])
end

C = [ [0.0, 1.0],[0.0, 0.0],[6.283185, 1.0],[6.283185, 0.0] ]
D = [1, 2, 1, 2]
Jbn = ARGV[0].to_i
Dp = 10

Bs = Bspline.new(XY + C, 7, 5, D);

vv = []
XY.each do |p|
	printf "% f, % f\n", p[0], p[1]
	vv.push p[0]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = Bs.plot(vv, 10, Jbn) do |a, b|
	printf "% f % f\n", a, b[0]
end
# Draw Graph
require "gnuplot"

Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'sin(x)'
		plot.ylabel 'Y'
		plot.xlabel 'X'
		x = vv.map {|v| v[0]}
		y = vv.map {|v| v[1]}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end
		y = x.map {|v| Bs.sekibun(v) - 1}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.title = "Integral"
		end
	end
end
