#! /usr/local/bin/ruby
require 'pp'
require '../pspline.so'
include PSPLINE

=begin

1 varable M dimension parametric interpolation

【Module】		PSPLINE
【Module Function】

（1）obj = PSPLINE#fft_complex_transform([[x1,y1],...,[xn,yn]], f)
	obj => FFT:[[u1,...,un],[v1,...,vn]]
	:1 list of complex data
	:2 type 1 => forward, 0 => backword, -1 => inverse

（2）obj = PSPLINE#fft_complex_get([[u1,...,un],[v1,...,vn]], t)
	obj => [x,y]
	:1 FFT_complex data
	:2 data point ( -N/2 <= t <= N/2 + N%2 )

（3）obj = PSPLINE#fft_complex_bspline([[u1,...,un],[v1,...,vn]], j)
	obj => PSPLINE#Bspline
	:1 FFT_complex data
	:2 dimension

=end

puts "# complex FFT data points"

XY = [ [0, 1],[0, 1],[0, 1],[1, 0,],[1, 0,],[1, 0],[1, 0],[0, 1],[0, 1],[0, 1] ]

fs = fft_complex_transform(XY, 1);
fs[0].each {|x| printf("% .2f ", x) }
puts
fs[1].each {|y| printf("% .2f ", y) }
puts

-5.upto(5) {|t|
	w = fft_complex_get(fs, t)
	printf "% .2f ", w[0]
}
puts
-5.upto(5) {|t|
	w = fft_complex_get(fs, t)
	printf "% .2f ", w[1]
}
puts

vv = []
bs = fft_complex_bspline(fs, 1) {|t, z| vv.push t}

puts "# Interpolation points"

s = bs.plot(vv, 4) do |a, b|
	printf "% .2f %10.7f %10.7f\n", a, b[0], b[1]
end
#STDERR.puts s

st = fft_complex_transform(fs, -1)
printf "["
st.each {|x| printf "% .1f ", x[0]}
printf "]\n["
st.each {|x| printf "% .1f ", x[1]}
puts "]"

