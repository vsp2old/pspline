#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

N varable M dimension parametric periodic interpolation

【Module】		PSPLINE
【Class】 		Pspline
【Method】
（1）new			Initialize
	obj = Pspline.new([[[y0,...],...],...[[yn-1,...],...]],[[x0,...,xn],...],[n,...],[j,...],[1,...])
	:1 list of data points.
	:2 list of variables
	:3 size of variables
	:4 dimension
	:5 type
（2）[]			Calculate interpolation
	obj[x,...]
（3）value		Calculate interpolation with differential value
	obj.value([x,...], b = 0)
	b:order of differential value （optional）
（4）plot
	obj.plot([[x0,...,xn],...], d, b = 0) { |[x,...],[y,...]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

# Parametric period interpolation

nst = [8, 3]
jst = [5, 2]
sst = [1, 0]

puts "# ２変数パラメトリックスプラインの補間計算 （閉曲面補間）"

S = [0, 1, 2, 3, 4, 5, 6, 7, 8]
T =	[0, 1, 2]

XYZ = [ [[ 1.000000, 0.000000, 0.0],[ 1.000000, 0.000000, 1.0],[ 1.000000, 0.000000, 2.0]],
		[[ 0.707107, 0.707107, 0.0],[ 0.707107, 0.707107, 1.0],[ 0.707107, 0.707107, 2.0]],
		[[ 0.000000, 1.000000, 0.0],[ 0.000000, 1.000000, 1.0],[ 0.000000, 1.000000, 2.0]],
		[[-0.707107, 0.707107, 0.0],[-0.707107, 0.707107, 1.0],[-0.707107, 0.707107, 2.0]],
		[[-1.000000, 0.000000, 0.0],[-1.000000, 0.000000, 1.0],[-1.000000, 0.000000, 2.0]],
		[[-0.707107,-0.707107, 0.0],[-0.707107,-0.707107, 1.0],[-0.707107,-0.707107, 2.0]],
		[[ 0.000000,-1.000000, 0.0],[ 0.000000,-1.000000, 1.0],[ 0.000000,-1.000000, 2.0]],
		[[ 0.707107,-0.707107, 0.0],[ 0.707107,-0.707107, 1.0],[ 0.707107,-0.707107, 2.0]] ]

Jbn = ARGV[0].to_i

Dp = 8

Qs = Pspline.new(XYZ, [S, T], nst, jst, sst)

for i in 0...8
	for j in 0..2
		printf("<%d,%d>", S[i], T[j])
		u = XYZ[i][j]
		printf(" % .7f % .7f % .7f\n", u[0], u[1], u[2])
	end
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
vv = [S, T]
N = Qs.plot(vv, Dp, [Jbn,0]) do |a, b|
	if (a[1] == 0.125)
		printf("%.3f %.3f % .7f % .7f % .7f\n", a[0], a[1], b[0], b[1], b[2])
	end
end
#STDERR.puts N
