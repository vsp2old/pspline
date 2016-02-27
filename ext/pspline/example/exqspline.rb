#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# parametric period interpolation (Closed curve interpolation)"

X = [0, 1, 2, 3, 4, 5, 6, 7, 8]

Y = [ [ 1.0, 0.0],[ 0.707107, 0.707107],[ 0.0, 1.0],[-0.707107, 0.707107],
	  [-1.0, 0.0],[-0.707107,-0.707107],[ 0.0,-1.0],[ 0.707107,-0.707107] ]

Jbn = ARGV[0].to_i

Dp = 8

Qs = Pspline.new(Y, [X], [8], [5], [1])

Y.each do |p|
	printf "% f, % f\n", p[0], p[1]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
vv = [X]
s = Qs.plot(vv, Dp, [Jbn]) do |a, b|
	printf "% .3f % f % f\n", a[0], b[0], b[1]
end
#STDERR.puts s
