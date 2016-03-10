/* Extention for ruby */

#include "ruby.h"

#ifdef __cplusplus
# define VALUEFUNC(f) ((VALUE (*)(...))f)
# define VOIDFUNC(f) ((void (*)(...))f)
#else
# define VALUEFUNC(f) (f)
# define VOIDFUNC(f) (f)
#endif

static VALUE
_wrap_ary_flatten(VALUE arg)
{
	VALUE result = rb_ary_new();
	if (TYPE(arg) == T_ARRAY) {
		if(RARRAY_LEN(arg) == 0)
			return result;
		else {
			VALUE last = rb_ary_pop(arg);
			return rb_ary_concat(_wrap_ary_flatten(arg), _wrap_ary_flatten(last));
		}
	} else 
		return rb_ary_push(result, arg);
}

#include "bspline.h"

static VALUE mPspline;

static VALUE cBspline;

static void Wrap_free_bspline(bspline<double> *arg0) { delete arg0; }

#define Wrap_bspline(klass, ptr) (\
	(ptr) ? Data_Wrap_Struct(klass, 0, VOIDFUNC(Wrap_free_bspline), ptr) : Qnil )

#define Get_bspline(val, ptr) {\
	if (NIL_P(val)) ptr = NULL;\
	else {\
	    if (!rb_obj_is_kind_of(val, cBspline))\
	        rb_raise(rb_eTypeError, "wrong argument type (expected bspline)");\
	    Data_Get_Struct(val, bspline<double>, ptr);\
	    if (!ptr) rb_raise(rb_eRuntimeError, "This bspline already released");\
	}\
}

//  bspline#new()

static VALUE
_wrap_new_bspline(int argc, VALUE *argv, VALUE klass)
{
	VALUE self = Data_Wrap_Struct(klass, 0, VOIDFUNC(Wrap_free_bspline), 0);
	rb_obj_call_init(self, argc, argv);
	return self;
}
/*******************************************************************************
	Calculate interpolation of 1 variable M dimension B-spline.
	bspline.new([...,[Xi,Yi0,...,YiM-1],...], N, J, S)
	  X0,..,Xi,..,XN-1:点のX座標
	  Yi0, ..., YiM-1:点のY座標
	  J:スプライン次数（階数K=J+1）
	Result:B-spline object
*******************************************************************************/

//  bspline#initialize

static VALUE
_wrap_init_bspline(int argc, VALUE *argv, VALUE self)
{
	VALUE va, varg, vargn, vargj, vargs, vargd = Qnil;
	int argp, argn, argj, args = 0, k = 1;
	int *argd = NULL;
	double *x, *y;

	rb_scan_args(argc, argv, "31", &varg, &vargn, &vargj, &vargs);

	Check_Type(varg, T_ARRAY);
	argp = RARRAY_LEN(varg);
	argn = NUM2INT(vargn);
	argj = NUM2INT(vargj);
	int D = argp - argn;
	if (D < 0)
		rb_raise(rb_eArgError, "Data points %d < Data size %d", argp, argn);
	if (D > argj)
		rb_raise(rb_eArgError, "Additional points %d > Degree %d", D, argj);
	if (argn <= argj)
		rb_raise(rb_eArgError, "Data points %d <= Degree %d", argn, argj);
	if (argc > 3) {
		vargd = argv[3];
		if (TYPE(vargd) == T_ARRAY) {
			int d = RARRAY_LEN(vargd);
			if (d != D)
				rb_raise(rb_eArgError, "Differential orders = %d, it should be %d", d, D);
			argd = ALLOC_N(int, argn + D);
			for (int i = 0; i < argn; ++i)
				argd[i] = 0;
			for (int i = 0; i < D; ++i)
				argd[argn + i] = NUM2INT(RARRAY_PTR(vargd)[i]);
		} else
			args = NUM2INT(vargs);
	}
	va = RARRAY_PTR(varg)[0];
	Check_Type(va, T_ARRAY);
	k = RARRAY_LEN(va) - 1;
	x = ALLOC_N(double, argp);
	y = ALLOC_N(double, argp * k);
	for (int i = 0; i < argp; i++) {
		if (i > 0) {
			va = RARRAY_PTR(varg)[i];
			Check_Type(va, T_ARRAY);
			if (RARRAY_LEN(va) - 1 != k)
				rb_raise(rb_eArgError, "Dimension should be %d", k);
		}
		for (int j = 0; j <= k; ++j)
			if (j == 0) x[i] = NUM2DBL(RARRAY_PTR(va)[0]);
			else  y[i*k+j-1] = NUM2DBL(RARRAY_PTR(va)[j]);
	}
try {
		Int s(Colon, argp - args%2, k);
		poly_array<double> argy(y, s);
//    	printf("n = %d, j = %d, s = %d\n", argn, argj, args);
//    	argy.print("\n");
		DATA_PTR(self) = new bspline<double>(argy, argn + args%2, x, argj, args, argd);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	free(y);
	free(x);
	if (argd) free(argd);
	return self;
}

static VALUE
make_list(int argc, VALUE *argv)
{
	VALUE vargs = rb_ary_new();
	for (int i = 0; i < argc; i++) {
		if (TYPE(argv[i]) == T_ARRAY) {
			VALUE val = rb_obj_dup(argv[i]);
			vargs = rb_ary_concat(vargs, _wrap_ary_flatten(val));
		}
		else
			rb_ary_push(vargs, argv[i]);
	}
	return vargs;
}
/*******************************************************************************
	関数値

	bspline#[x[, ...]]
	  x:点のX座標（数値，または数値の配列）

	Result:補間値
*******************************************************************************/
static VALUE
_wrap_bspline_bracket(int argc, VALUE *argv, VALUE self)
{
	bspline<double> *bsp;
	double arg;
	VALUE val, vargs;
	VALUE vresult = Qnil;

	Get_bspline(self, bsp);
	int K = bsp->Unit_size();

	vargs = make_list(argc, argv);
	argc = RARRAY_LEN(vargs);
try {
	vresult = rb_ary_new();
	for (int i = 0; i < argc; i++) {
		arg = NUM2DBL(RARRAY_PTR(vargs)[i]);
		poly_array<double> result = (*bsp)[arg];
		if (K == 1)
			val = rb_float_new(result[0]);
		else {
			val = rb_ary_new();
			for (int j = 0; j < K; ++j)
				rb_ary_push(val, rb_float_new(result[j]));
		}
		rb_ary_push(vresult, val);
	}
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	return argc == 1 ? rb_ary_shift(vresult) : vresult;
}
/*******************************************************************************
	関数値
	bspline#value(x, d = 0)
	  x:点のX座標
	  d:微分階数（省略可）
	Result:補間値，または微分値
*******************************************************************************/
static VALUE
_wrap_bspline_value(int argc, VALUE *argv, VALUE self)
{
	VALUE varg1, varg2 ;
	bspline<double> *arg0 ;
	double arg1 ;
	int arg2 = 0 ;
	VALUE vresult = Qnil;
	
	rb_scan_args(argc, argv, "11", &varg1, &varg2);
	Get_bspline(self, arg0);
	int K = arg0->Unit_size();
	arg1 = NUM2DBL(varg1);
	if (argc > 1) arg2 = NUM2INT(varg2);
try {
	poly_array<double> result = (*arg0)(arg1, arg2);
	if (K == 1) 
		vresult = rb_float_new(result[0]);
	else {
		vresult = rb_ary_new();
		for (int j = 0; j < K; ++j)
			rb_ary_push(vresult, rb_float_new(result[j]));
	}
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	return vresult;
}

static VALUE
_wrap_bspline_sekibun(VALUE self, VALUE varg)
{
	bspline<double> *arg0 ;
	VALUE vresult = Qnil;

	Get_bspline(self, arg0);
	int K = arg0->Unit_size();
try {
	poly_array<double> result = arg0->sekibun(NUM2DBL(varg));
	if (K == 1) 
		vresult = rb_float_new(result[0]);
	else {
		vresult = rb_ary_new();
		for (int j = 0; j < K; ++j)
			rb_ary_push(vresult, rb_float_new(result[j]));
	}
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	return vresult;
}

static VALUE
_wrap_bspline_plot(int argc, VALUE *argv, VALUE self)
{
	bspline<double> *arg0;		// bspline object
	VALUE varg1;				// list of point
	VALUE varg3; int arg3;		// number of division
	VALUE varg4; int arg4 = 0;	// order differential
	rb_scan_args(argc, argv, "21", &varg1, &varg3, &varg4);
	Check_Type(varg1, T_ARRAY);
	arg3 = NUM2INT(varg3);
	if (argc > 2) arg4 = NUM2INT(varg4);
	Get_bspline(self, arg0);
	int argn = RARRAY_LEN(varg1);
	double *argx = ALLOC_N(double, argn);
	for (int i = 0; i < argn; ++i) {
		VALUE x = rb_ary_shift(varg1);
		argx[i] = NUM2DBL(x);
	}
	int K = arg0->Unit_size();
	double *arg1, *arg2[K];
	int result;
try {
	result = plot<double>(*arg0, argn, argx, &arg1, arg2, arg3, arg4);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	for (int i = 0; i < result; i++) {
		VALUE v1 = rb_float_new(arg1[i]);
		VALUE v2 = rb_ary_new();
		for (int k = 0; k < K; ++k)
			rb_ary_push(v2, rb_float_new(arg2[k][i]));
		VALUE vv = rb_ary_new();
		rb_ary_push(vv, v1);
		rb_ary_push(vv, v2);
		rb_ary_push(varg1, vv);
		if (rb_block_given_p()) rb_yield(vv);
	}
	free(arg1);
	for (int i = 0; i < K; ++i) free(arg2[i]);
	free(argx);
	return INT2NUM(result);
}

static VALUE
_wrap_complex_trans(VALUE self, VALUE varg1, VALUE varg2)
{
	VALUE varg, vargx, vargy, vresult = Qnil;
	Check_Type(varg1, T_ARRAY);
	int argn = RARRAY_LEN(varg1);
	int argf = NUM2INT(varg2);
	if (argf != 1) {
		if (argn != 2) rb_raise(rb_eArgError, "Argument is not Complex");
		vargx = RARRAY_PTR(varg1)[0]; Check_Type(vargx, T_ARRAY);
		vargy = RARRAY_PTR(varg1)[1]; Check_Type(vargy, T_ARRAY);
		 argn = RARRAY_LEN(vargx);
		if (argn != RARRAY_LEN(vargy)) rb_raise(rb_eArgError, "Argument is not Complex");
	}
	double *dfft = ALLOC_N(double, argn << 1);
	double *data = ALLOC_N(double, argn << 1);
	if (argf == 1) {
		for (int i = 0; i < argn; i++) {
			varg = RARRAY_PTR(varg1)[i];
			Check_Type(varg, T_ARRAY);
			if (RARRAY_LEN(varg) != 2) rb_raise(rb_eArgError, "Argument is not Complex");
			dfft[2*i] = NUM2DBL(RARRAY_PTR(varg)[0]);
			dfft[2*i+1] = NUM2DBL(RARRAY_PTR(varg)[1]);
		}
try {
	fft<>::complex_transform(dfft, argn, argf);
	fft<>::half_en_pack(data, argn, dfft);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
		vresult = rb_ary_new();
		vargx = rb_ary_new();
		vargy = rb_ary_new();
		rb_ary_push(vresult, vargx);
		rb_ary_push(vresult, vargy);
		for (int i = 0; i < argn; i++) {
			rb_ary_push(vargx, rb_float_new(data[i]));
			rb_ary_push(vargy, rb_float_new(data[i+argn]));
		}
	} else {
		for (int i = 0; i < argn; i++) {
			data[i] = NUM2DBL(RARRAY_PTR(vargx)[i]);
			data[i+argn] = NUM2DBL(RARRAY_PTR(vargy)[i]);
		}
try {
	fft<>::half_de_pack(dfft, argn, data);
	fft<>::complex_transform(dfft, argn, argf);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
		vresult = rb_ary_new();
		for (int i = 0; i < argn; i++) {
			varg = rb_ary_new();
			rb_ary_push(vresult, varg);
			rb_ary_push(varg, rb_float_new(dfft[2*i]));
			rb_ary_push(varg, rb_float_new(dfft[2*i+1]));
		}
	}
	free(data);
	free(dfft);
	return vresult;
}

static VALUE
_wrap_complex_get(VALUE self, VALUE varg1, VALUE varg2)
{
	VALUE vargx, vargy, vresult = Qnil;
	Check_Type(varg1, T_ARRAY);
	int argn = RARRAY_LEN(varg1);
	int args = NUM2INT(varg2);
	if (argn != 2) rb_raise(rb_eArgError, "Argument is not Complex");
	vargx = RARRAY_PTR(varg1)[0];
	vargy = RARRAY_PTR(varg1)[1];
	 argn = RARRAY_LEN(vargx);
	if (RARRAY_LEN(vargy) != argn) rb_raise(rb_eArgError, "Argument is not Complex");
	double *data = ALLOC_N(double, argn << 1);
	for (int i = 0; i < argn; i++) {
		data[i] = NUM2DBL(RARRAY_PTR(vargx)[i]);
		data[i+argn] =NUM2DBL(RARRAY_PTR(vargy)[i]);
	}
	double w[2];
	fft<>::half_get(data, argn, args, w);
	vresult = rb_ary_new();
	rb_ary_push(vresult, rb_float_new(w[0]));
	rb_ary_push(vresult, rb_float_new(w[1]));
	free(data);
	return vresult;
}

static VALUE
_wrap_real_trans(VALUE self, VALUE varg1, VALUE varg2)
{
	VALUE vresult = Qnil;
	Check_Type(varg1, T_ARRAY);
	int argn = RARRAY_LEN(varg1);
	int argf = NUM2INT(varg2);
	double *data = ALLOC_N(double, argn);
	for (int i = 0; i < argn; i++)
		data[i] = NUM2DBL(RARRAY_PTR(varg1)[i]);
try {
	fft<>::real_transform(data, argn, argf);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	vresult = rb_ary_new();
	for (int i = 0; i < argn; i++)
		rb_ary_push(vresult, rb_float_new(data[i]));
	free(data);
	return vresult;
}

static VALUE
_wrap_real_get(VALUE self, VALUE varg1, VALUE varg2)
{
	VALUE vresult = Qnil;
	Check_Type(varg1, T_ARRAY);
	int argn = RARRAY_LEN(varg1);
	int args = NUM2INT(varg2);
	double *data = ALLOC_N(double, argn);
	for (int i = 0; i < argn; i++)
		data[i] = NUM2DBL(RARRAY_PTR(varg1)[i]);
	double w[2];
	fft<>::real_get(data, argn, args, w);
	vresult = rb_ary_new();
	rb_ary_push(vresult, rb_float_new(w[0]));
	rb_ary_push(vresult, rb_float_new(w[1]));
	free(data);
	return vresult;
}

static VALUE
_wrap_complex_bspline(VALUE self, VALUE varg1, VALUE varg2)
{
	Check_Type(varg1, T_ARRAY);
	int argn = RARRAY_LEN(varg1);
	if (argn != 2) rb_raise(rb_eArgError, "Argument is not Complex");
	VALUE vargx = RARRAY_PTR(varg1)[0];
	VALUE vargy = RARRAY_PTR(varg1)[1];
	Check_Type(vargx, T_ARRAY);
	Check_Type(vargy, T_ARRAY);
	argn = RARRAY_LEN(vargx);
	if (RARRAY_LEN(vargy) != argn) rb_raise(rb_eArgError, "Argument is not Complex");
	double *data = ALLOC_N(double, argn << 1);
	for (int i = 0; i < argn; i++) {
		data[i] = NUM2DBL(RARRAY_PTR(vargx)[i]);
		data[i+argn] = NUM2DBL(RARRAY_PTR(vargy)[i]);
	}
	double *dfft = ALLOC_N(double, argn << 1);
	fft<>::half_de_pack(dfft, argn, data);
	int argj = NUM2INT(varg2);
	bspline<> *result;
try {
	result = fft<>::cfft(dfft, argn, argj);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	if (rb_block_given_p()) {
		double val[2]; int N = argn / 2;
		for (int i = -N; i <= N + argn%2; i++) {
			fft<>::half_get(data, argn, i, val);
			VALUE vv = rb_ary_new();
			rb_ary_push(vv, INT2NUM(i));
			VALUE xy = rb_ary_new();
			rb_ary_push(xy, rb_float_new(val[0]));
			rb_ary_push(xy, rb_float_new(val[1]));
			rb_ary_push(vv, xy);
			rb_yield(vv);
		}
	}
	free(dfft);
	free(data);
	return Wrap_bspline(cBspline, result);
}

static VALUE
_wrap_real_bspline(VALUE self, VALUE varg1, VALUE varg2)
{
	Check_Type(varg1, T_ARRAY);
	int argn = RARRAY_LEN(varg1);
	double *data = ALLOC_N(double, argn);
	for (int i = 0; i < argn; i++)
		data[i] = NUM2DBL(RARRAY_PTR(varg1)[i]);
	int argj = NUM2INT(varg2);
	bspline<> *result;
try {
	result = fft<>::rfft(data, argn, argj);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	if (rb_block_given_p()) {
		double val[2]; int N = argn / 2;
		for (int i = -N; i <= N + argn%2; i++) {
			fft<>::real_get(data, argn, i, val);
			VALUE vv = rb_ary_new();
			rb_ary_push(vv, INT2NUM(i));
			VALUE xy = rb_ary_new();
			rb_ary_push(vv, xy);
			rb_ary_push(xy, rb_float_new(val[0]));
			rb_ary_push(xy, rb_float_new(val[1]));
			rb_yield(vv);
		}
	}
	free(data);
	return Wrap_bspline(cBspline, result);
}

static VALUE vfunc;

static double ufunc(poly_array<double>& pa)
{
	int K = pa.unit_size();
	VALUE *vp = ALLOC_N(VALUE, K);
	for (int i = 0; i < K; i++) vp[i] = rb_float_new(pa[i]);
	VALUE vresult = rb_funcall2(vfunc, rb_intern("call"), K, vp);
	free(vp);
	return NUM2DBL(vresult);
}

static VALUE
_wrap_bspline_lint(int argc, VALUE *argv, VALUE self)
{
	bspline<double> *arg0;		// bspline object
	VALUE vargu, vargx, vargy;
	Get_bspline(self, arg0);
	rb_scan_args(argc, argv, "22", &vfunc, &vargu, &vargx, &vargy);
	int c = arg0->Icox(), j = arg0->Jisu(), n;
	double *q = arg0->Knots(), x = arg0->x_min(), y = arg0->x_max();
	Check_Type(vargu, T_ARRAY);
	int argn = RARRAY_LEN(vargu);
	n = argn > 0 ? argn : c - j + 1;
	double *u = ALLOC_N(double, n);
	for (int i = 0; i < n; i++)
		if (argn > 0) {
			u[i] = NUM2DBL(RARRAY_PTR(vargu)[i]);
			if (u[i] < x || y < u[i] || (i > 0 && u[i] <= u[i-1]))
				rb_raise(rb_eArgError, "Illegal argument");
		} else {
			u[i] = x + i * (y - x) / (n - 1);
			rb_ary_push(vargu, rb_float_new(u[i]));
		}
	x = (argc >= 3) ? NUM2DBL(vargx) : u[0];
	y = (argc == 4) ? NUM2DBL(vargy) :   0 ;
	bspline<double> *result;
try {
	result = arg0->line_integral(ufunc, n, u, x, y);
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	free(u);
	return Wrap_bspline(cBspline, result);
}

static VALUE _wrap_bspline_mul(VALUE self, VALUE varg);

static void
_wrap_Init_bspline(void)
{
	cBspline = rb_define_class_under(mPspline, "Bspline", rb_cObject);
	rb_include_module(cBspline, rb_mEnumerable);
	rb_define_singleton_method(cBspline, "new", VALUEFUNC(_wrap_new_bspline), -1);
	rb_define_method(cBspline, "initialize", VALUEFUNC(_wrap_init_bspline), -1);
	rb_define_method(cBspline, "[]", VALUEFUNC(_wrap_bspline_bracket), -1);
	rb_define_method(cBspline, "value", VALUEFUNC(_wrap_bspline_value), -1);
	rb_define_method(cBspline, "sekibun", VALUEFUNC(_wrap_bspline_sekibun), 1);
	rb_define_method(cBspline, "plot", VALUEFUNC(_wrap_bspline_plot), -1);
	rb_define_method(cBspline, "line_integral", VALUEFUNC(_wrap_bspline_lint), -1);
	rb_define_method(cBspline, "*", VALUEFUNC(_wrap_bspline_mul), 1);
	rb_define_module_function(mPspline, "fft_complex_transform", VALUEFUNC(_wrap_complex_trans), 2);
	rb_define_module_function(mPspline, "fft_complex_get", VALUEFUNC(_wrap_complex_get), 2);
	rb_define_module_function(mPspline, "fft_complex_bspline", VALUEFUNC(_wrap_complex_bspline), 2);
	rb_define_module_function(mPspline, "fft_real_transform", VALUEFUNC(_wrap_real_trans), 2);
	rb_define_module_function(mPspline, "fft_real_get", VALUEFUNC(_wrap_real_get), 2);
	rb_define_module_function(mPspline, "fft_real_bspline", VALUEFUNC(_wrap_real_bspline), 2);
}

static VALUE cPspline;

static void Wrap_free_pspline(pspline<double> *arg0) { delete arg0; }

#define Wrap_pspline(klass, ptr) (\
	(ptr) ? Data_Wrap_Struct(klass, 0, VOIDFUNC(Wrap_free_pspline), ptr) : Qnil )

#define Get_pspline(val, ptr) {\
	if (NIL_P(val)) ptr = NULL;\
	else {\
	    if (!rb_obj_is_kind_of(val, cPspline))\
	        rb_raise(rb_eTypeError, "wrong argument type (expected pspline)");\
	    Data_Get_Struct(val, pspline<double>, ptr);\
	    if (!ptr) rb_raise(rb_eRuntimeError, "This pspline already released");\
	}\
}

//	Pspline#new()

static VALUE
_wrap_new_pspline(int argc, VALUE *argv, VALUE klass)
{
	VALUE self = Data_Wrap_Struct(klass, 0, VOIDFUNC(Wrap_free_pspline), 0);
	rb_obj_call_init(self, argc, argv);
	return self;
}

#define Get_param(N, varg, ptr, str) {\
	Check_Type(varg, T_ARRAY);\
	int k = RARRAY_LEN(varg);\
	if (k != N)\
		rb_raise(rb_eArgError, "%s list length = %d, it should be %d", str, k, N);\
	ptr = ALLOC_N(int, k);\
	for (int i = 0; i < k; ++i) ptr[i] = NUM2INT(RARRAY_PTR(varg)[i]);\
}

static comb_array<double> make_var(int N, VALUE var, int *n = NULL, int *m = NULL)
{
	int M = 0, s[N]; VALUE vs[N];
	for (int i = 0; i < N; ++i) {
		vs[i] = RARRAY_PTR(var)[i];
		Check_Type(vs[i], T_ARRAY);
		s[i] = RARRAY_LEN(vs[i]);
		if (n != NULL && m != NULL) {
			int vars = n[i] + m[i]%2;
			if (s[i] != vars)
				rb_raise(rb_eArgError, "Elements = %d, it should be %d", s[i], vars);
		}
		M += s[i];
	}
	double v[M], *p = v;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < s[i]; ++j)
			*p++ = NUM2DBL(RARRAY_PTR(vs[i])[j]); 
	return comb_array<double>(v, N, s, 1);
}

static VALUE make_flat(int K, VALUE val, int N, int *n)
{
	if (N == 0) return val;
	if (K != n[0]) rb_raise(rb_eArgError, "Wrong argument Array");
	VALUE vals = rb_ary_new();
	for (int i = 0; i < K; ++i) {
		VALUE va = RARRAY_PTR(val)[i];
		if (TYPE(va) == T_ARRAY) {
			VALUE vb = rb_obj_dup(va);
			vals = rb_ary_concat(vals, make_flat(RARRAY_LEN(va), vb, N-1, n+1));
		}
		else
			rb_ary_push(vals, va);
	}
	return vals;
}

static poly_array<double> make_val(int &K, VALUE val, int N, int *n)
{
	if (K != n[0]) rb_raise(rb_eArgError, "Member size = %d, it should be %d", K, n[0]);
	VALUE vals = rb_ary_new();
	for (int i = 0; i < K; ++i) {
		VALUE va = RARRAY_PTR(val)[i];
		if (TYPE(va) == T_ARRAY) {
			VALUE vb = rb_obj_dup(va);
			vals = rb_ary_concat(vals, make_flat(RARRAY_LEN(va), vb, N-1, n+1));
		}
		else
			rb_ary_push(vals, va);
	}
	int M = RARRAY_LEN(vals);
	K = M;
	for (int i = 0; i < N; ++i) K /= n[i];
	double v[M];
	for (int i = 0; i < M; ++i) v[i] = NUM2DBL(RARRAY_PTR(vals)[i]);
	return poly_array<double>(v, N, n, K);
}
/*******************************************************************************
	Calculate interpolation of N variable M dimension P-spline.
	Pspline.new([X^Y^...], [S,T,...], N, J, S)
	  X^Y^...:点の座標
	  S,T,...:変数の座標
	  N:
	  J:スプライン次数（階数K=J+1）
	  S:
	Result:P-spline object
*******************************************************************************/

//	Pspline#initialize

static VALUE
_wrap_init_pspline(int argc, VALUE *argv, VALUE self)
{
	VALUE vxyz, vstu, vargn, vargj, vargs;

	int N, K, *n, *j, *s;

	rb_scan_args(argc, argv, "5", &vxyz, &vstu, &vargn, &vargj, &vargs);

	Check_Type(vstu, T_ARRAY);
	N = RARRAY_LEN(vstu);
	
	Get_param(N, vargs, s, "Type");

	Get_param(N, vargj, j, "Degree");

	Get_param(N, vargn, n, "Size");

	Check_Type(vxyz, T_ARRAY);
	K = RARRAY_LEN(vxyz);

try {

	comb_array<double> stu = make_var(N, vstu, n, s);
	poly_array<double> xyz = make_val(K, vxyz, N, n);

	Int argn(N, n, K);
	Int argj(N, j, 0);
	Int args(N, s, 0);
//		argn.print(", "); argj.print(", "); args.print("\n");
//		stu.print("\n"); xyz.print("\n");
	DATA_PTR(self) = new pspline<double>(xyz, stu, argj, args);

} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	free(n);
	free(j);
	free(s);
	return self;
}
/*******************************************************************************
	関数値

	pspline#[x[, ...]]
	  x:点のX座標（数値，または数値の配列）

	Result:補間値
*******************************************************************************/
static VALUE
_wrap_pspline_bracket(int argc, VALUE *argv, VALUE self)
{
	pspline<double> *bsp;
	VALUE vargs, vresult, mresult = Qnil;

	Get_pspline(self, bsp);
	int N = bsp->Grid();
	int K = bsp->Unit_size();
	mresult = rb_ary_new();
try {
	for (int m = 0; m < argc; ++m) {
		vargs = argv[m];
		poly_array<double> val;
		if (TYPE(vargs) == T_ARRAY) {
			int L = RARRAY_LEN(vargs);
			if (L != N)
				rb_raise(rb_eArgError, "Array length = %d, it should be %d.", L, N);
			double args[N+1];
			for (int i = 0; i < N; ++i)
				args[i] = NUM2DBL(RARRAY_PTR(vargs)[i]);
			poly<double> arg(args, N);
			val = (*bsp)[arg];
		} else
			val = (*bsp)[NUM2DBL(vargs)];
		if (K == 1) 
			vresult = rb_float_new(val[0]);
		else {
			vresult = rb_ary_new();
			for (int j = 0; j < K; ++j)
				rb_ary_push(vresult, rb_float_new(val[j]));
		}
		rb_ary_push(mresult, vresult);
	}
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	return argc == 1 ? rb_ary_shift(mresult) : mresult;
}

void ps_eArgError(VALUE err)
{
	VALUE str = rb_String(err);
	rb_raise(rb_eArgError, "wrong argument %s", StringValuePtr(str));
}
/*******************************************************************************
	関数値
	pspline#value(x, d = 0)
	  x:点のX座標
	  d:微分階数（省略可）
	Result:補間値，または微分値
*******************************************************************************/
static VALUE
_wrap_pspline_value(int argc, VALUE *argv, VALUE self)
{
	VALUE varg1, varg2, varg3;
	pspline<double> *arg0;
	int argn, arg2 = 0, *carg2 = NULL;
	double arg1, *args, *carg3 = NULL;
	VALUE vargs, vresult = Qnil;

	rb_scan_args(argc, argv, "12", &varg1, &varg2, &varg3);

	Get_pspline(self, arg0);
	int N = arg0->Grid();
	int K = arg0->Unit_size();
	// 全微分
	if (argc == 3) {
		arg2 = NUM2INT(varg2);
		carg3 = new double[N];
		Check_Type(varg3, T_ARRAY);
		if (RARRAY_LEN(varg3) != N) ps_eArgError(varg3);
		for (int i = 0; i < N; ++i)
			carg3[i] = NUM2DBL(RARRAY_PTR(varg3)[i]);
	}
	// 偏微分
	if (argc == 2) {
		if (TYPE(varg2) == T_ARRAY && RARRAY_LEN(varg2) == N) {
			carg2 = new int[N];
			for (int i = 0; i < N; ++i)
				carg2[i] = NUM2INT(RARRAY_PTR(varg2)[i]);
		} else if (TYPE(varg1) != T_ARRAY) {
			arg2 = NUM2INT(varg2);
			argc = 0;
		} else ps_eArgError(varg2);
	}
	// Argument list
	poly<double> arg;
	args = NULL;
	if (TYPE(varg1) == T_ARRAY) {
		argn = RARRAY_LEN(varg1);
		vargs = make_list(argn, RARRAY_PTR(varg1));
		argn = RARRAY_LEN(vargs);
		if (argn == N) {
			args = new double[N+1];
			for (int i = 0; i < N; ++i)
				args[i] = NUM2DBL(RARRAY_PTR(vargs)[i]);
		} else ps_eArgError(varg1);
	} else
		if (argc <= 1) arg1 = NUM2DBL(varg1);
		else ps_eArgError(varg1);
	if (argc < 2) argc ^= 1;
	if (args) arg = poly<double>(args, N);
	poly_array<double> val;
try {
	switch (argc) {
		case 0: 	// 関数値
			val = args ? (*arg0)(arg) : (*arg0)(arg1);
			delete[] args;
			break;
		case 1: 	// 関数積の微分
			val = (*arg0)(arg1, arg2);
			delete[] args;
			break;
		case 2: 	// 偏微分 (partial derivative)
			val = (*arg0)(arg, carg2);
			delete[] carg2;
			break;
		case 3: 	// 全微分 (total derivative)
			val = (*arg0)(arg, arg2, carg3);
			delete[] carg3;
	}
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	if (K == 1) 
		vresult = rb_float_new(val[0]);
	else {
		vresult = rb_ary_new();
		for (int j = 0; j < K; ++j)
			rb_ary_push(vresult, rb_float_new(val[j]));
	}
	return vresult;
}

static VALUE
_wrap_pspline_sekibun(int argc, VALUE *argv, VALUE self)
{
	pspline<double> *bsp;
	VALUE vargs, vresult = Qnil;

	Get_pspline(self, bsp);
	int N = bsp->Grid();
	int K = bsp->Unit_size();

	vargs = make_list(argc, argv);
	argc = RARRAY_LEN(vargs);
try {
	if (argc == N || argc == 1) {
		double s[N+1];
		for (int i = 0; i < N; ++i)
			if (i == 0 || argc == N)
				s[i] = NUM2DBL(RARRAY_PTR(vargs)[i]);
			else
				s[i] = s[i-1];
		poly<double> args(s, N);
		poly_array<double> val = bsp->sekibun(args);
		if (K == 1) 
			vresult = rb_float_new(val[0]);
		else {
			vresult = rb_ary_new();
			for (int j = 0; j < K; ++j)
				rb_ary_push(vresult, rb_float_new(val[j]));
		}
	} else rb_raise(rb_eArgError, "wrong argument Sekibun");
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	return vresult;
}

static VALUE
_wrap_pspline_plot(int argc, VALUE *argv, VALUE self)
{
	pspline<double> *arg0;
	VALUE vstu, varg, vdp, vjbn;
	int arg3, *arg4 = NULL;

	rb_scan_args(argc, argv, "21", &varg, &vdp, &vjbn);

	Get_pspline(self, arg0);
	int N = arg0->Grid();
	int K = arg0->Unit_size();
	double *arg1[N], *arg2[K];
	arg3 = NUM2INT(vdp);
	if (argc > 2) {
		arg4 = new int[N];
		if (TYPE(vjbn) == T_ARRAY) {
			int n = RARRAY_LEN(vjbn);
			if (n != N)
				rb_raise(rb_eArgError, "Differential orders = %d, it should be %d", n, N);
			for (int i = 0; i < N; ++i)
				arg4[i] = NUM2INT(RARRAY_PTR(vjbn)[i]);
		} else
			arg4[0] = NUM2INT(vjbn);
	}
	Check_Type(varg, T_ARRAY);
	N = RARRAY_LEN(varg);
	vstu = rb_ary_new();
	int result;
	for (int i = 0; i < N; ++i) rb_ary_push(vstu, rb_ary_shift(varg));
try {
	if (TYPE(RARRAY_PTR(vstu)[0]) == T_ARRAY) {
		comb_array<double> stu = make_var(N, vstu);
		result = plot<double>(*arg0, stu, arg1, arg2, arg3, arg4);
	} else {
		double *x = new double[N];
		for (int i = 0; i < N; i++) x[i] = NUM2DBL(RARRAY_PTR(vstu)[i]);
		result = plot<double>(*arg0, N, x, arg1, arg2, arg3, arg4[0]);
		delete[] x;
		N = 1;
	}
} catch (const char *c) {
	rb_raise(rb_eRuntimeError, "%s", c);
}
	delete[] arg4;
	for (int i = 0; i < result; i++) {
		VALUE v1 = rb_ary_new();
		for (int j = 0; j < N; ++j) rb_ary_push(v1, rb_float_new(arg1[j][i]));
		VALUE v2 = rb_ary_new();
		for (int k = 0; k < K; ++k) rb_ary_push(v2, rb_float_new(arg2[k][i]));
		VALUE vv = rb_ary_new();
		rb_ary_push(vv, v1);
		rb_ary_push(vv, v2);
		rb_ary_push(varg, vv);
		if (rb_block_given_p()) rb_yield(vv);
	}
	for (int i = 0; i < N; ++i) FREE(arg1[i]);
	for (int i = 0; i < K; ++i) FREE(arg2[i]);
	return INT2NUM(result);
}

static VALUE
_wrap_bspline_mul(VALUE self, VALUE varg)
{
	bspline<double> *arg0;
	Get_bspline(self, arg0);
	pspline<double> *result;
	if (rb_obj_is_kind_of(varg, cBspline)) {
		bspline<double> *arg1;
	    Data_Get_Struct(varg, bspline<double>, arg1);
	    if (!arg1) rb_raise(rb_eArgError, "Not Bspline object");
		result = new pspline<double>((*arg0) * (*arg1));
	} else if (rb_obj_is_kind_of(varg, cPspline)) {
		pspline<double> *arg1;
	    Data_Get_Struct(varg, pspline<double>, arg1);
	    if (!arg1) rb_raise(rb_eArgError, "Not Pspline Object");
		result = new pspline<double>((*arg0) * (*arg1));
	} else {
		double a = NUM2DBL(varg);
		bspline<double> *result = new bspline<double>(*arg0);
		double *p = *result;
		int k = result->Unit_size();
		int c = ((base_spline<double>*)result)->Icox();
		for (int i = 0; i < c * k; i++) p[i] *= a;
		return Wrap_bspline(cBspline, result);
	}
	return Wrap_pspline(cPspline, result);
}

static VALUE
_wrap_pspline_mul(VALUE self, VALUE varg)
{
	pspline<double> *arg0;
	Get_pspline(self, arg0);
	pspline<double> *result;
	if (rb_obj_is_kind_of(varg, cBspline)) {
		bspline<double> *arg1;
	    Data_Get_Struct(varg, bspline<double>, arg1);
	    if (!arg1) rb_raise(rb_eArgError, "Not Bspline object");
		result = new pspline<double>((*arg0) * (*arg1));
	}/* else if (rb_obj_is_kind_of(varg, cPspline)) {
		pspline<double> *arg1;
	    Data_Get_Struct(varg, pspline<double>, arg1);
	    if (!arg1) rb_raise(rb_eArgError, "Not Pspline object");
		result = new pspline<double>((*arg0) * (*arg1));
	}*/ else {
		double a = NUM2DBL(varg);
		result = new pspline<double>(*arg0);
		double *p = *result;
		const Int& c = result->values->icox;
		int s = c.grid_size();
		for (int i = 0; i < s; i++) p[i] *= a;
	}
	return Wrap_pspline(cPspline, result);
}

static void
_wrap_Init_pspline(void)
{
	cPspline = rb_define_class_under(mPspline, "Pspline", rb_cObject);
	rb_include_module(cPspline, rb_mEnumerable);
	rb_define_singleton_method(cPspline, "new", VALUEFUNC(_wrap_new_pspline), -1);
	rb_define_method(cPspline, "initialize", VALUEFUNC(_wrap_init_pspline), -1);
	rb_define_method(cPspline, "[]", VALUEFUNC(_wrap_pspline_bracket), -1);
	rb_define_method(cPspline, "value", VALUEFUNC(_wrap_pspline_value), -1);
	rb_define_method(cPspline, "sekibun", VALUEFUNC(_wrap_pspline_sekibun), -1);
	rb_define_method(cPspline, "plot", VALUEFUNC(_wrap_pspline_plot), -1);
	rb_define_method(cPspline, "*", VALUEFUNC(_wrap_pspline_mul), 1);
}

#ifdef __cplusplus
extern "C" {
#endif

void Init_pspline(void)
{
	mPspline = rb_define_module("PSPLINE");

	_wrap_Init_bspline();
	_wrap_Init_pspline();
}

#ifdef __cplusplus
}
#endif

