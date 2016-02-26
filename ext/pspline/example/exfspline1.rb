#! /usr/local/bin/ruby
require '../pspline.so'
include PSPLINE

=begin

1 varable M dimension parametric interpolation

【Module】		PSPLINE
【Class】 		Bspline
【Method】
（1）new		Initialize
	obj = Pspline.new([[x1,y1,...],...,[xn,yn,...]], n, j)
	:1 list of data points.
	:2 size of data points
	:3 dimension
（2）[]			Calculate interpolation
	obj[x]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b: order of differential value （optional）
（4）plot
	obj.plot([x1,...,xn], d, b = 0) { |x,[y,...]| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

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

