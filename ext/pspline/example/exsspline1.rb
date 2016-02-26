#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

1 varable M dimension Riesenfeld periodic interpolation

【Module】		PSPLINE
【Class】 		Bspline
【Method】
（1）new		Initialize
	obj = Bspline.new([[x0, y0,...],...,[xn, yn,...]], n, j, 3)
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

#puts "# Riesenfeld interpolation （Closed curve interpolation）"

n = 8
j = 3
s = 3

# リーゼンフェルトスプラインの補間計算（閉曲線補間）

XY = [ [0, 4.0, 0.0],[1, 3.0, 3.0],[2, 0.0, 4.0],[3,-3.0, 3.0],
	   [4,-4.0, 0.0],[5,-3.0,-3.0],[6, 0.0,-4.0],[7, 3.0,-3.0],
	   [8, 4.0, 0.0] ]

Jbn = ARGV[0].to_i

Dp = 10

Sp = Bspline.new(XY, n, j, s)

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
s = Sp.plot(vv, Dp, Jbn) do |a, b|
	printf "% .2f % .7f % .7f\n", a, b[0], b[1]
end
#STDERR.puts s
