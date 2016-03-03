#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# Parametric period interpolation (Closed curve)"

XY = [ [ 0, 1.0, 0.0],[ 1, 0.923880, 0.382683],[ 2, 0.707107, 0.707107],[ 3, 0.382683, 0.923880],
	   [ 4, 0.0, 1.0],[ 5,-0.382683, 0.923880],[ 6,-0.707107, 0.707107],[ 7,-0.923880, 0.382683],
	   [ 8,-1.0, 0.0],[ 9,-0.923880,-0.382683],[10,-0.707107,-0.707107],[11,-0.382683,-0.923880],
	   [12, 0.0,-1.0],[13, 0.382683,-0.923880],[14, 0.707107,-0.707107],[15, 0.923880,-0.382683],
	   [16, 1.0, 0.0] ]

Jbn = ARGV[0].to_i

Dp = 8

Qs = Bspline.new(XY, 16, 9, 1)

vv = []
XY.each do |p|
	printf "% f, % f\n", p[1], p[2]
	vv.push p[0]
end

Sint = Qs.line_integral(lambda{|x,y|Math.sqrt(x*x+y*y)*0.5}, vv)
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = Qs.plot(vv, Dp, Jbn) do |a, b|
	printf "%6.3f % f % f % f\n", a, b[0], b[1], Sint[a]
end
#STDERR.puts s
