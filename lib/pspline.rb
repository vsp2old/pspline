require "pspline/version"
require "pspline.so"

module PSPLINE

class Rfft

	def initialize(data)
		@R = fft_real_transform(data, 1)
	end
	def [] (t)
		fft_real_get(@R, t)
	end
	def spline(j)
		@V = []
		fft_real_bspline(@R, j) {|t, z| @V.push t}
	end
	def axis
		@V.dup
	end
	def inverse
		fft_real_transform(@R, -1)
	end
	def backword
		fft_real_transform(@R, 0)
	end
	def real
		@R.dup
	end

end

class Cfft

	def initialize(data)
		@F = fft_complex_transform(data, 1)
	end
	def [] (t)
		fft_complex_get(@F, t)
	end
	def spline(j)
		@V = []
		fft_complex_bspline(@F, j) {|t, z| @V.push t}
	end
	def inverse
		fft_complex_transform(@F, -1)
	end
	def backword
		fft_real_transform(@F, 0)
	end
	def axis
		@V.dup
	end
	def real
		@F[0].dup
	end
	def imag
		@F[1].dup
	end

end

end
