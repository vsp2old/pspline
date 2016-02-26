#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

1 varable M dimension parametric interpolation

【Module】		PSPLINE
【Class】 		Pspline
【Method】
（1）new		Initialize
	obj = Pspline.new([[y1,...],...,[yn,...]],[[x1,...,xn]],[n],[j],[0])
	:1 list of data points.
	:2 list of variables
	:3 size of variables
	:4 dimension
	:5 type
（2）[]			Calculate interpolation
	obj[x]
（3）value		Calculate interpolation with differential value
	obj.value([x], [b])
	b: order of differential value （optional）
	obj.value([x], b, [d])
	d: differential vector
（4）plot
	obj.plot([[x1,...,xn]], d, [b] = nil) { |[x],[y,...]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

puts "# Parametric interpolation"

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
