#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

puts "# complex FFT data points"

XY = [ [0, 1],[0, 1],[0, 1],[1, 0,],[1, 0,],[1, 0],[1, 0],[0, 1],[0, 1],[0, 1] ]

fs = Cfft.new(XY)
fs.real.each {|x| printf("% .2f ", x) }
puts
fs.image.each {|y| printf("% .2f ", y) }
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

bs = fs.spline order:1
vv = fs.axis

puts "# Interpolation points"

s = bs.plot(vv, 4) do |a, b|
	printf "% .2f %10.7f %10.7f\n", a, b[0], b[1]
end
#STDERR.puts s

st = fs.inverse
printf "["
st.each {|x| printf "% .1f ", x[0]}
printf "]\n["
st.each {|x| printf "% .1f ", x[1]}
puts "]"

