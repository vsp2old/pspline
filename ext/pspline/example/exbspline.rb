#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

Ad = [[-4.0, 1.34095e-3],[-3.6, 2.98189e-3],[-3.2, 6.62420e-3],[-2.8, 1.46827e-2],
	  [-2.4, 3.23838e-2],[-2.0, 7.06508e-2],[-1.6, 1.50527e-1],[-1.2, 3.05020e-1],
	  [-0.8, 5.59055e-1],[-0.4, 8.55639e-1],[ 0.4, 8.55639e-1],[ 0.8, 5.59055e-1],
	  [ 1.2, 3.05020e-1],[ 1.6, 1.50527e-1],[ 2.0, 7.06508e-2],[ 2.4, 3.23838e-2],
	  [ 2.8, 1.46827e-2],[ 3.2, 6.62420e-3],[ 3.6, 2.98189e-3],[ 4.0, 1.34095e-3]]

Jbn = ARGV[0].to_i

Dp = 10

Bs = Bspline.new(Ad, 20, 5)

puts "# Interpolation of the sec2(x) function"

vv = []
Ad.each do |p|
	printf "% .2f % f\n", p[0], p[1]
	vv.push p[0];
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
elsif Jbn > 0
	printf ", Jbn = %d\n", Jbn
elsif Jbn < 0
	printf ", Jsk = %d\n", -Jbn
end

s = Bs.plot(vv, Dp, Jbn) do |a,b|
	printf "% .2f % f\n", a, b[0]
end
# STDERR.puts s
