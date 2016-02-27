#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

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
	printf "%7.4f, % f\n", p[0], p[1]
	vv.push p[0]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = Bs.plot(vv, 10, Jbn) do |u, v|
	printf "%8.5f % f\n", u, v[0]
end
