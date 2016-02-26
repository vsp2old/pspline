#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

N variable M dimension Reisenfeld interpolation

【Module】		PSPLINE
【Class】 		Pspline
【Method】
（1）new		Initialize
	obj = Pspline.new([[[y1,...],...],...,[[yn,...],...]],[[x1,...,xn],...],[n,...],[j,...],[2,...])
	:1 list of data points.
	:2 list of variables
	:3 size of variables
	:4 dimension
	:5 type
（2）[]			Calculate interpolation
	obj[x,...]
（3）value		Calculate interpolation with differential value
	:1 obj.value([x,...], b = 0)
	:2 order of differential value （optional）
（4）plot
	obj.plot([[x1,...,xn],...], d, [b,...]) { |[x,...],[y,...]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

# Reisenfeld interpolation

nn = [4, 6]
jj = [3, 2]
ss = [3, 2]

puts "# リーゼンフェルトスプライン曲面の補間計算"

S = [0, 1, 2, 3, 4]
T = [0, 1, 2, 3, 4, 5]

XYZ = [ [[ 0.0, 5.0, 9.0],[ 0.0, 5.0, 6.0],[ 0.0, 0.5, 6.0],[ 0.0, 0.5, 2.0],[ 0.0, 3.0, 1.0],[ 0.0, 3.0, 0.0]],
		[[ 5.0,-2.5, 6.5],[ 5.0,-2.5, 3.5],[ 0.5,-0.25,5.75],[0.5,-0.25,1.75],[3.0,-1.5,-0.5],[ 3.0,-1.5,-1.5]],
		[[ 0.0,-5.0, 9.0],[ 0.0,-5.0, 6.0],[ 0.0,-0.5, 6.0],[ 0.0,-0.5, 2.0],[ 0.0,-3.0, 1.0],[ 0.0,-3.0, 0.0]],
		[[-5.0, 2.5,11.5],[-5.0, 2.5, 8.5],[-0.5,0.25,6.25],[-0.5,0.25,2.25],[-3.0, 1.5, 2.5],[-3.0, 1.5, 1.5]] ]

Jbn = ARGV[0].to_i

Dp = 8

Rs = Pspline.new(XYZ, [S, T], nn, jj, ss)

for i in 0..3
	for j in 0..5
		printf("<%d,%d>", S[i], T[j])
		u = XYZ[i][j]
		printf(" % .7f % .7f %10.7f\n", u[0], u[1], u[2])
	end
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
vv = [S, T]
N = Rs.plot(vv, Dp, [Jbn,0]) do |a, b|
	if (a[1] == 0.0)
		printf("%.3f %.3f % .7f % .7f %10.7f\n", a[0], a[1], b[0], b[1], b[2])
	end
end
#STDERR.puts N
