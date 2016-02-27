#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# Parametric spline interpolation"

X = [0,1,2,3,4,5,6,7]
Y = [ [ 1.0, 0.0], [ 0.0, 1.0], [-1.0, 0.0], [ 0.0,-1.0],
	  [ 2.0, 0.0], [ 0.0, 2.0], [-2.0, 0.0], [ 0.0,-2.0] ]
Jbn = ARGV[0].to_i
Dp = 20

Ps = Pspline.new(Y, [X], [8], [2], [0])

Y.each do |p|
	printf "% .2f % .2f\n", p[0], p[1]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
vv = [X]
s = Ps.plot(vv, Dp, [Jbn]) do |a, b|
	c = Ps.sekibun(a[0])
	printf "% .2f % .7f % .7f % .7f % .7f\n", a[0], b[0], b[1], c[0], c[1]
end
#STDERR.puts s
