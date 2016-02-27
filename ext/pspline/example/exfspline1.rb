#! /usr/local/bin/ruby
require 'pspline'
include PSPLINE

printf "# real FFT data points\n"

X = [ 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 ]

fs = fft_real_transform(X, 1);
fs.each {|x| printf("% .2f ", x) }
puts

printf "Real :"
-5.upto(5) {|t|
	w = fft_real_get(fs, t)
	printf "% .2f ", w[0]
}
puts
printf "Imag :"
-5.upto(5) {|t|
	w = fft_real_get(fs, t)
	printf "% .2f ", w[1]
}
puts

vv = []
bs = fft_real_bspline(fs, 3) {|t, z| vv.push t}

puts "# Interpolation points"

s = bs.plot(vv, 4) do |a, b|
	printf "% .2f %10.7f %10.7f\n", a, b[0], b[1]
end
#STDERR.puts s

st = fft_real_transform(fs, -1)
printf "["
st.each {|x| printf "% .1f ", x}
puts "]"

