#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

1 varable M dimension parametric interpolation

【Module】		PSPLINE
【Class】 		Bspline
【Method】
（1）new		Initialize
	obj = Bspline.new([[x1,y1,...],...,[xn,yn,...]], n, j)
	:1 list of data points.
	:2 size of data points
	:3 dimension
（2）[]			Calculate interpolation
	obj[x]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b: order of differential value （optional）
（4）plot
	obj.plot([x1,...,xn], d, b = 0) { |x,[y,...]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

puts "# Parametric interpolation"

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
