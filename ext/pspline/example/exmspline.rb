#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

Ad = [[-4.0, 1.34095e-3],[-3.6, 2.98189e-3],[-3.2, 6.62420e-3],[-2.8, 1.46827e-2],
	  [-2.4, 3.23838e-2],[-2.0, 7.06508e-2],[-1.6, 1.50527e-1],[-1.2, 3.05020e-1],
	  [-0.8, 5.59055e-1],[-0.4, 8.55639e-1],[ 0.4, 8.55639e-1],[ 0.8, 5.59055e-1],
	  [ 1.2, 3.05020e-1],[ 1.6, 1.50527e-1],[ 2.0, 7.06508e-2],[ 2.4, 3.23838e-2],
	  [ 2.8, 1.46827e-2],[ 3.2, 6.62420e-3],[ 3.6, 2.98189e-3],[ 4.0, 1.34095e-3]]

Bd = [[-4.0, -4.0],[-3.6, -3.6],[-3.2, -3.2],[-2.8, -2.8],
	  [-2.4, -2.4],[-2.0, -2.0],[-1.6, -1.6],[-1.2, -1.2],
	  [-0.8, -0.8],[-0.4, -0.4],[ 0.4,  0.4],[ 0.8,  0.8],
	  [ 1.2,  1.2],[ 1.6,  1.6],[ 2.0,  2.0],[ 2.4,  2.4],
	  [ 2.8,  2.8],[ 3.2,  3.2],[ 3.6,  3.6],[ 4.0,  4.0]]

Jbn = ARGV[0].to_i

Dp = 10

As = Bspline.new(Ad, 20, 5)
Bs = Bspline.new(Bd, 20, 2)
Ps = As * Bs

puts "# Interpolation of multiple Bspline function"

vv = []
Ad.each do |p|
	printf "% .2f % f\n", p[0], p[1]
	vv.push p[0];
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
elsif Jbn > 0
	printf ", Jbn = %d\n", Jbn
elsif Jbn < 0
	printf ", Jsk = %d\n", -Jbn
end

s = Ps.plot(vv, Dp, Jbn) do |a,b|
	c = As.sekibun(a[0])
	d = Ps.sekibun(a[0])
	printf "% .2f % f % f % f\n", a[0], b[0], c, d
end
# STDERR.puts s
require "gnuplot"
 
Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'Multiple Bspline'
		plot.ylabel 'Y'
		plot.xlabel 'X'
 		x = vv.map {|u,v| u[0]}
		y = vv.map {|u,v| v[0]}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end
		y = vv.map {|u,v| Ps[[u[0], 0.5*u[0]]]}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end
	end
end
