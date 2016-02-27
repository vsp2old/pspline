#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# Parametric period interpolation (Closed curve interpolation)"

XY = [ [0, 1.0, 0.0],[1, 0.707107, 0.707107],[2, 0.0, 1.0],[3,-0.707107, 0.707107],
	   [4,-1.0, 0.0],[5,-0.707107,-0.707107],[6, 0.0,-1.0],[7, 0.707107,-0.707107],
	   [8, 1.0, 0.0] ]

Jbn = ARGV[0].to_i

Dp = 8

Qs = Bspline.new(XY, 8, 5, 1)

vv = []
XY.each do |p|
	printf "% f, % f\n", p[1], p[2]
	vv.push p[0]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = Qs.plot(vv, Dp, Jbn) do |a, b|
	printf "% .2f % f % f\n", a, b[0], b[1]
end
#STDERR.puts s
