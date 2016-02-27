#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# Riesenfeld interpolation （Closed curve interpolation）"

n = 8
j = 3
s = 3

XY = [ [0, 4.0, 0.0],[1, 3.0, 3.0],[2, 0.0, 4.0],[3,-3.0, 3.0],
	   [4,-4.0, 0.0],[5,-3.0,-3.0],[6, 0.0,-4.0],[7, 3.0,-3.0],
	   [8, 4.0, 0.0] ]

Jbn = ARGV[0].to_i

Dp = 10

Sp = Bspline.new(XY, n, j, s)

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
s = Sp.plot(vv, Dp, Jbn) do |a, b|
	printf "% .2f % .7f % .7f\n", a, b[0], b[1]
end
#STDERR.puts s
