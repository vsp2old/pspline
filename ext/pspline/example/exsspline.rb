#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

n = 8
j = 3
s = 3

puts "# Riesenfeld interpolation （Closed curve）"

X = [0, 1, 2, 3, 4, 5, 6, 7, 8]

Y = [ [ 4.0, 0.0],[ 3.0, 3.0],[ 0.0, 4.0],[-3.0, 3.0],
	  [-4.0, 0.0],[-3.0,-3.0],[ 0.0,-4.0],[ 3.0,-3.0] ]

Jbn = ARGV[0].to_i

Dp = 10

Sp = Pspline.new(Y, [X], [n], [j], [s])

Y.each do |p|
	printf "% f, % f\n", p[0], p[1]	
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
vv = [X]
s = Sp.plot(vv, Dp, [Jbn]) do |a, b|
	printf "% .2f, <% .7f, % .7f>\n", a[0], b[0], b[1]
end
#STDERR.puts s
