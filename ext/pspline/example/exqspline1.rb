#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

1 varable M dimension parametric periodic interpolation

【Module】		PSPLINE
【Class】 		Bspline
【Method】
（1）new			Initialize
	obj = Bspline.new([[x0,y0,...],...,[xn,yn,...]], n, j, 1)
	:1 list of data points.
	:2 size of data points
	:3 dimension
	:4 type
（2）[]			Calculate interpolation
	obj[x]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b:order of differential value （optional）
（4）plot
	obj.plot([x0,...,xn], d, b = 0) { |x,[y,...]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

#puts "# Parametric interpolation （Closed surface interpolation）"

XY = [ [0, 1.0, 0.0],[1, 0.707107, 0.707107],[2, 0.0, 1.0],[3,-0.707107, 0.707107],
	   [4,-1.0, 0.0],[5,-0.707107,-0.707107],[6, 0.0,-1.0],[7, 0.707107,-0.707107],
	   [8, 1.0, 0.0] ]

Jbn = ARGV[0].to_i

Dp = 8

Qs = Bspline.new(XY, 8, 5, 1)

vv = []
XY.each do |p|
#	printf "%f, %f\n", p[1], p[2]
	vv.push p[0]
end
#printf "# value of interpolation points, Dp = %d", Dp
#if Jbn == 0
#	print "\n"
#else
#	printf ", Jbn = %d\n", Jbn
#end
s = Qs.plot(vv, Dp, Jbn) do |a, b|
	printf "% .2f % .15f % .15f\n", a, b[0], b[1]
end
#STDERR.puts s
