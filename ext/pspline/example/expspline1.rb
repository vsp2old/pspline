#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# Parametric spline interpolation"

XY = [ [0, 1.0, 0.0], [1, 0.0, 1.0], [2,-1.0, 0.0], [3, 0.0,-1.0],
	   [4, 2.0, 0.0], [5, 0.0, 2.0], [6,-2.0, 0.0], [7, 0.0,-2.0] ]
Jbn = ARGV[0].to_i
Dp = 20

Qs = Bspline.new(XY, 8, 2)

vv = []
XY.each do |p|
	printf "<% .1f % .1f>\n", p[1], p[2]
	vv.push p[0]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = Qs.plot(vv, Dp, Jbn) do |a, b|
	c = Qs.sekibun(a);
	printf "% .2f % .7f % .7f % .7f % .7f\n", a, b[0], b[1], c[0], c[1]
end
#STDERR.puts s
