#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

=begin

Interpolation with boundary condition by additional data points

【Module】		PSPLINE
【Class】 		Bspline
【Method】
（1）new			Initialize
	obj = Bspline.new([[x1,y1],...,[xn,yn],[xn+1,yn+1],...,[xn+d,yn+d]]], n, j)
	:1 list of data points.
	:2 number of data points
	:3 dimension
（2）[]			Calculate interpolation
	obj[x]
	obj[x0,...,xi]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b:order of differential value （optional）
（4）plot
	obj.plot([x0,...,xn], d, b = 0) { |x,[y]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

n = 7
j = 5

#
#puts "# Interpolation of the Bessel function"
#
#	Data points
Ad = [ [0.0, 1.0    ],[0.8, 0.84629],[1.6, 0.45540],[2.0, 0.22389],
	   [2.4, 0.00251],[3.2,-0.32019],[4.0,-0.39715], ]
#	Additional data points
C = [ [0.4, 0.96040],[1.2, 0.67113],[2.8,-0.18504],[3.6,-0.39177] ]
Jbn = ARGV[0].to_i
Dp = 10

Bs = Bspline.new(Ad + C, n, j)

vv = []
for i in 0..6
	(p, q) = Ad[i]
#	printf "%.2f, %f\n", p, q
	vv.push p
end
#printf "# value of interpolation points, Dp = %d", Dp
#if Jbn == 0
#	print "\n"
#else
#	printf ", Jbn = %d\n", Jbn
#end
s = Bs.plot(vv, Dp, Jbn) do |u, v|
	printf "% .15f % .15f\n", u, v[0]
end
#puts s
