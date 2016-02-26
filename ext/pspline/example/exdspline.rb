#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

Interpolation with boundary condition by differential value

【Module】		PSPLINE
【Class】 		Bspline
【Method】
（1）new		Initialize
	obj = Bspline.new([[x1,y1],...,[xn,yn],[xn+1,yn+1],...,[xn+d,yn+d]],n,j,[b1,...,bd])
	:1 list of data points.
	:2 number of data points
	:3 dimension
	:4 list of order of differential value
（2）[]			Calculate interpolation
	obj[x]
	obj[x0,...,xi]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b: order of differential value （optional）
（4）plot
	obj.plot([x0,...,xn], d, b = 0) { |x,[y]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

puts "# Interpolation of the sin(x) function"

XY = [ [0.0],[1.0],[2.0],[3.0],[4.0],[5.0],[6.283185] ]

XY.each do |a|
	a.push Math.sin(a[0])
end

C = [ [0.0, 1.0],[0.0, 0.0],[6.283185, 1.0],[6.283185, 0.0] ]
D = [1, 2, 1, 2]
Jbn = ARGV[0].to_i
Dp = 10

Bs = Bspline.new(XY + C, 7, 5, D);

vv = []
XY.each do |p|
	printf "% f, % f\n", p[0], p[1]
	vv.push p[0]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = Bs.plot(vv, 10, Jbn) do |a, b|
	printf "% .7f % .15f\n", a, b[0]
end
