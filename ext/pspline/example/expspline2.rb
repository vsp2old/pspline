#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

# Parametric interpolation

nn = [5, 3]
jj = [3, 2]
ss = [0, 0]

puts "# 2 variable 3 dimension parametric spline interpolation"

S = [0, 1, 2, 3, 4]
T = [0, 1, 2]

XYZ = [ [[ 2.0,-1.0,-2.0], [ 2.0,-2.0,-1.0], [ 2.0,-2.0,-3.0]],
		[[ 1.0, 0.0,-1.0], [ 1.0,-1.0, 0.0], [ 1.0,-1.0,-2.0]],
		[[ 0.0, 2.0, 0.0], [ 0.0, 0.0, 2.0], [ 0.0, 0.0,-2.0]],
		[[-1.0, 2.0, 1.0], [-1.0, 1.0, 2.0], [-1.0, 1.0, 0.0]],
		[[-2.0, 3.0, 2.0], [-2.0, 2.0, 3.0], [-2.0, 2.0, 1.0]] ]

Jbn = ARGV[0].to_i

Dp = 10

Ps = Pspline.new(XYZ, [S, T], nn, jj, ss)

for i in 0..4
	for j in 0..2
		printf "%d,%d", S[i], T[j]
		u = XYZ[i][j]
		printf " <% .1f % .1f % .1f>\n", u[0], u[1], u[2]
	end
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
vv = [S, T]
N = Ps.plot(vv, Dp, [Jbn,0]) do |a, b|
	if (a[0] == 0.1)
		printf "% .2f % .2f % f % f % f\n", a[0], a[1], b[0], b[1], b[2]
	end
end
#STDERR.puts N
