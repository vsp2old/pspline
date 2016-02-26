require 'mkmf'

dir_config('PSPLINE','.')
$CPPFLAGS += " -std=c++11"
if have_header('bspline.h')
	create_makefile('pspline')
end
