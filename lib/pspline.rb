require "pspline/version"
require "pspline.so"

module PSPLINE

class Rfft < Array

	def initialize(data)
		super(fft_real_transform(data, 1))
	end
	def [] (t)
		fft_real_get(self, t)
	end
	def spline order:
		@V = []
		fft_real_bspline(self, order) {|t, z| @V.push t}
	end
	def axis
		@V.dup
	end
	def inverse
		fft_real_transform(self, -1)
	end
	def backword
		fft_real_transform(self, 0)
	end

end

class Cfft < Array

	def initialize(data)
		super(fft_complex_transform(data, 1))
	end
	def [] (t)
		fft_complex_get(self, t)
	end
	def spline order:
		@V = []
		fft_complex_bspline(self, order) {|t, z| @V.push t}
	end
	def inverse
		fft_complex_transform(self, -1)
	end
	def backword
		fft_real_transform(self, 0)
	end
	def axis
		@V.dup
	end
	def real
		self.first.dup
	end
	def image
		self.last.dup
	end

end

end

class Array

	def fft_complex_forward
		PSPLINE::Cfft.new(self)
	end
	def fft_real_forward
		PSPLINE::Rfft.new(self)
	end

end
