#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

1 varable M dimension Riesenfeld interpolation

【Module】		PSPLINE
【Class】 		Pspline
【Method】
（1）new		Initialize
	obj = Pspline.new([[y1,...],...,[yn,...]],[[x1,...,xn]],[n],[j],[2])
	:1 list of data points.
	:2 list of variables
	:3 size of variables
	:4 dimension
	:5 type
（2）[]			Calculate interpolation
	obj[x]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b:order of differential value （optional）
（4）plot
	obj.plot([[x1,...,xn]], d, b = 0) { |[x],[y,...]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

n   = 7
j   = 2
s   = 2

puts "# Riesenfeld interpolation"

X = [0, 1, 2, 3, 4, 5, 6]
Y = [ [ 0.0, 0.0],[ 3.0, 0.0],[ 3.0, 1.0],[ 0.5, 2.0],
	  [ 0.5, 6.0],[ 5.0, 6.0],[ 5.0, 9.0] ]

Jbn = ARGV[0].to_i

Dp = 10

Rs = Pspline.new(Y, [X], [n], [j], [s])
	
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
s = Rs.plot(vv, Dp, [Jbn]) do |a, b|
	printf "% .2f, <% .7f, % .7f>\n", a[0], b[0], b[1]
end
#STDERR.puts s
