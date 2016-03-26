#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

printf "# real FFT data points\n"

X = [ 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 ]

fs = Rfft.new(X);
fs.each {|x| printf("% .2f ", x) }
puts

printf " Real:"
-5.upto(5) {|t|
	w = fs[t]
	printf "% .2f ", w[0]
}
puts
printf "Image:"
-5.upto(5) {|t|
	w = fs[t]
	printf "% .2f ", w[1]
}
puts

bs = fs.spline order:3
vv = fs.axis

puts "# Interpolation points"

s = bs.plot(vv, 4) do |a, b|
	printf "% .2f %10.7f %10.7f\n", a, b[0], b[1]
end
#STDERR.puts s

st = fs.inverse
printf "["
st.each {|x| printf "% .1f ", x}
puts "]"

