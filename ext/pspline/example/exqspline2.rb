#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

# Parametric period interpolation

nst = [8, 3]
jst = [5, 2]
sst = [1, 0]

puts "# 2 variable parametric period interpolation (Closed surface interpolation)"

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
		printf(" % f % f % f\n", u[0], u[1], u[2])
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
		printf("%.3f %.3f % f % f % f\n", a[0], a[1], b[0], b[1], b[2])
	end
end
#STDERR.puts N
