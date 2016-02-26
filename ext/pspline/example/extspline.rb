#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

Interpolation with periodic boundary condition

【Module】		PSPLINE
【Class】 		Bspline
【Method】
（1）new			Initialize
	obj = Bspline.new([[x0,y0],...,[xn,yn]], n, j, 1)
	:1 list of data points
	:2 number of data points
	:3 dimension
	:4 flag indicating the periodic function
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

#
#puts "# Interpolation of the Jacobi function"
#
XY = [ [ 0.0,    0.0    ],[ 0.5474, 0.5    ],[ 1.2837, 0.86603],[ 2.7681, 1.0    ],
	   [ 4.2525, 0.86603],[ 4.9888, 0.5    ],[ 5.5362, 0.0    ],[ 6.0836,-0.5    ],
	   [ 6.8199,-0.86603],[ 8.3043,-1.0    ],[ 9.7887,-0.86603],[10.5250,-0.5    ],
	   [11.0724, 0.0    ] ]

Jbn = ARGV[0].to_i

Dp = 10

Bs = Bspline.new(XY, 12, 5, 1)

vv = []
XY.each do |p|
#	printf "%.5f, %f\n", p[0], p[1]
	vv.push p[0]
end
#printf "# value of interpolation points, Dp = %d", Dp
#if Jbn == 0
#	print "\n"
#else
#	printf ", Jbn = %d\n", Jbn
#end
s = Bs.plot(vv, 10, Jbn) do |u, v|
	printf "%10.7f % .7f\n", u, v[0]
end
