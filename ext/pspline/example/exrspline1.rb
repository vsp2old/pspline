#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

n   = 7
j   = 2
s   = 2

puts "# Riesenfeld interpolation"

XY = [ [0, 0.0, 0.0],[1, 3.0, 0.0],[2, 3.0, 1.0],[3, 0.5, 2.0],
	   [4, 0.5, 6.0],[5, 5.0, 6.0],[6, 5.0, 9.0] ]

Jbn = ARGV[0].to_i

Dp = 10

Rs = Bspline.new(XY, n, j, s)

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
s = Rs.plot(vv, Dp, Jbn) do |a, b|
	printf "% .2f % .7f % .7f\n", a, b[0], b[1]
end
#STDERR.puts s
