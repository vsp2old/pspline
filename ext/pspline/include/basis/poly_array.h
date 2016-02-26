#ifndef _POLY_ARRAY_5_H_INCLUDED_
#define _POLY_ARRAY_5_H_INCLUDED_
#include <stdarg.h>

template<typename T> class poly_array;
template<typename T> class poly_view;
template<typename T> class poly_array_lu;
template<typename T> class comb_array;
template<typename T> class comb_view;
template<typename T> class comp_array;
template<typename T> class comp_view;

typedef enum {Atom, Colon, Slate, Cube, Tetra, Penta, Hexa} grid_type;

template <typename T>
void poly_subst(int n, T *w)
{
	for (int i = 0; i <= n; ++i) w[i] = 1;
}

template <typename T, class... Args>
void poly_subst(int n, T *w, T car, Args... cdr)
{
	w[0] = car; if (n > 0) poly_subst(n-1, &w[1], cdr...);
}

/*******************************************************************************
	poly<T>(N, T[N]);
	[T0,...,TN]
*******************************************************************************/
template<typename T> class poly
{
template<class> friend class poly;
template<class> friend class poly_array;
template<class> friend class poly_view;
template<class> friend class poly_array_lu;
template<class> friend class comb_array;
template<class> friend class comb_view;
template<class> friend class comp_array;
template<class> friend class comp_view;

int N, *cp; T *w;	// data fields

public:
poly(): cp(NULL), w(NULL) {}
// constructor
//	poly<T>(T*, N);  alias of T*
poly(T *a, int n = 1) : N(n), cp(NULL), w(a) {} // [T0,...,Tn]
// poly<T>(N, T0, ...);
template <class... Args>
poly(int n, const T car, const Args... cdr) : N(n), cp(new int(1)), w(new T[n+1])
{
	w[0] = car; if (n > 0) poly_subst(n-1, &w[1], cdr...);
}
poly(int n, const T *s) : N(n), cp(new int(1)), w(new T[n+1])
{
	memcpy(w, s, (N + 1) * sizeof(T));
}
poly(int n, const T *s, T t) : N(n), cp(new int(1)), w(new T[n+1])
{
	memcpy(w, s, N * sizeof(T)); w[N] = t;
}
// copy construtor
template <class S>
poly(const poly<S>& p) : N(p.N), cp(new int(1)), w(new T[p.N+1])
{
	for (int i = 0; i <= N; ++i) w[i] = p.w[i];
}
poly(const poly& p) : N(p.N), cp(new int(1)), w(new T[p.N+1])
{
	for (int i = 0; i <= N; ++i) w[i] = p.w[i];
}
// destructor
~poly() { if (cp && --*cp == 0) { delete cp; delete[] w; }}
// substituton operator
poly& operator=(const poly& p)
{
	if (this != &p)
	{
		if (cp && --*cp == 0) { delete cp; delete[] w; }
		N = p.N; cp = new int(1); w = new T[p.N+1];
		for (int i = 0; i <= N; ++i) w[i] = p.w[i];
	}
	return *this;
}
// substitute grid size
poly& operator=(const T *s)
{
	memcpy(w, s, N * sizeof(T));
	return *this;
}
// substitute unit size
poly& operator=(T t)
{
	w[N] = t;
	return *this;
}
poly operator * (const poly& p) const
{
	assert(N == p.N);
	T result[N];
	for (int i=0;i<N;++i) result[i] = w[i] * p.w[i];
	return poly(N, result, w[N]);
}
poly operator / (const poly& p) const
{
	assert(N == p.N);
	T result[N];
	for (int i=0;i<N;++i) result[i] = w[i] / p.w[i];
	return poly(N, result, w[N]);
}
poly operator + (const poly& p) const
{
	assert(N == p.N);
	T result[N];
	for (int i=0;i<N;++i) result[i] = w[i] + p.w[i];
	return poly(N, result, w[N]);
}
poly operator - (const poly& p) const
{
	assert(N == p.N);
	T result[N];
	for (int i=0;i<N;++i) result[i] = w[i] - p.w[i];
	return poly(N, result, w[N]);
}
poly& operator *= (const poly& p)
{
	assert(N == p.N);
	for (int i=0;i<N;++i) w[i] *= p.w[i];
	return *this;
}
poly& operator /= (const poly& p)
{
	assert(N == p.N);
	for (int i=0;i<N;++i) w[i] /= p.w[i];
	return *this;
}
poly& operator += (const poly& p)
{
	assert(N == p.N);
	for (int i=0;i<N;++i) w[i] += p.w[i];
	return *this;
}
poly& operator -= (const poly& p)
{
	assert(N == p.N);
	for (int i=0;i<N;++i) w[i] -= p.w[i];
	return *this;
}
poly& operator %= (const poly& p)
{
	assert(N == p.N);
	for (int i=0;i<N;++i) w[i] %= p.w[i];
	return *this;
}
poly operator - () const
{
	T s[N];
	for (int i=0;i<N;++i) s[i] = - w[i];
	return poly(N, s, w[N]);
}
poly operator * (T a) const
{
	T s[N];
	for (int i=0;i<N;++i) s[i] = w[i] * a;
	return poly(N, s, w[N]);
}
poly operator / (T a) const
{
	T s[N];
	for (int i=0;i<N;++i) s[i] = w[i] / a;
	return poly(N, s, w[N]);
}
poly operator + (T a) const
{
	T s[N];
	for (int i=0;i<N;++i) s[i] = w[i] + a;
	return poly(N, s, w[N]);
}
poly operator - (T a) const
{
	T s[N];
	for (int i=0;i<N;++i) s[i] = w[i] - a;
	return poly(N, s, w[N]);
}
poly operator % (int a) const
{
	T s[N];
	for (int i=0;i<N;++i) s[i] = w[i] % a;
	return poly(N, s, w[N]);
}
poly& operator *= (T a)
{
	for(int i=0;i<N;i++) w[i] *= a;
	return *this;
}
poly& operator /= (T a)
{
	for(int i=0;i<N;i++) w[i] /= a;
	return *this;
}
poly& operator += (T a)
{
	for(int i=0;i<N;i++) w[i] += a;
	return *this;
}
poly& operator -= (T a)
{
	for(int i=0;i<N;i++) w[i] -= a;
	return *this;
}
poly& operator %= (T a)
{
	for(int i=0;i<N;i++) w[i] %= a;
	return *this;
}
bool operator == (const poly& p) const
{
	assert(N == p.N);
	int i; for(i=0;i<N;++i) if(w[i] != p.w[i]) break;
	return i == N && w[N] == p.w[N];
}
bool operator != (const poly& p) const
{
	assert(N == p.N);
	int i; for(i=0;i<N;i++) if(w[i] != p.w[i]) break;
	return i != N || w[N] != p.w[N];
}
bool operator > (const poly& p) const
{
	assert(N == p.N);
	int i; for(i=0;i<N;++i) if(w[i] <= p.w[i]) break;
	return i == N && w[N] == p.w[N];
}
bool operator >= (const poly& p) const
{
	assert(N == p.N);
	int i; for(i=0;i<N;++i) if(w[i] < p.w[i]) break;
	return i == N && w[N] == p.w[N];
}
bool operator < (const poly& p) const
{
	assert(N == p.N);
	int i; for(i=0;i<N;i++) if(w[i] >= p.w[i]) break;
	return i == N && w[N] == p.w[N];
}
bool operator <= (const poly& p) const
{
	assert(N == p.N);
	int i; for(i=0;i<N;i++) if(w[i] > p.w[i]) break;
	return i == N && w[N] == p.w[N];
}
bool operator == (int a) const
{
	int i;
	for (i = 0; i < N; ++i) if (w[i] != a) break;
	return i == N;
}
bool operator != (int a) const
{
	int i;
	for (i = 0; i < N; ++i) if (w[i] != a) break;
	return i != N;
}
bool operator < (int a) const
{
	int i;
	for (i = 0; i < N; ++i) if (w[i] >= a) break;
	return i == N;
}
bool operator >= (int a) const
{
	int i;
	for (i = 0; i < N; ++i) if (w[i] < a) break;
	return i == N;
}
bool operator > (int a) const
{
	int i;
	for (i = 0; i < N; ++i) if (w[i] <= a) break;
	return i == N;
}
bool operator <= (int a) const
{
	int i;
	for (i = 0; i < N; ++i) if (w[i] > a) break;
	return i == N;
}
// member function
int length() const {return N+1;}
// Π w[i] : i = n,...,N
T data_size(int n = 0) const
{
	assert(n <= N);
	T z = 1;
	for(int i = N; i > n; --i) z *= w[i-1];
	return z * w[N];
}
// (Σ com[i] : i = n,...,N-1) * unit_size 
T comb_size(int n = 0) const
{
	assert(n <= N);
	int z = (n == N ? 1 : 0);
	for(int i = N; i > n; --i) z += w[i-1];
	return z * w[N];
}
T cast_size(int n, int m = 0) const
{
	assert(0 < N && n <= N && m <= n);
	int z = 0;
	for(int i = m; i < n; ++i) z += w[i];
	return z * w[N];
}
// Π w[i] : i = 0,...,N-1
T grid_size() const
{
	T z = 1;
	for (int k = 0; k < N; ++k) z *= w[k];
	return z;
}
// Σ w[i] : i = 0,...,N-1
T band_size() const
{
	T z = 0;
	for (int k = 0; k < N; ++k) z += w[k];
	return z;
}
int grid() const {return N;}
//   w[N]
T unit_size() const { return w[N]; }
// pointer of w
operator T* () const { return w; }
// calculate Atom Position
T operator && (const poly& p) const
{
	int M = p.N;
	assert(N >= M);
	T sum = p.w[M];
	T sek = w[N];
	for (int i = N-1; i >= 0; --i) {
		if (i < M)
		sum += sek * p.w[i];
		sek *= w[i];
	}
	return sum;
}
bool is_variable() const
{
	int N = this->grid(); T n = this->data_size();
	return N == 0 && n == 1;
}
// dimension of grid
bool is_atom()  const { return N == Atom;  }
bool is_colon() const { return N == Colon; }
bool is_slate() const { return N == Slate; }
bool is_cube()  const { return N == Cube;  }
// type conversion operator
// Referrece of element.
T operator[] (int i) const { return w[i]; }
// operators
poly operator<<(int i) const { return poly(N-i, &w[i]); }
poly operator>>(int i) const { return poly(N-i, w, data_size(N-i)); }
poly operator+() { return poly(Colon, grid_size(), unit_size()); }
char *to_s() const
{
	char buf[512];
	char *c = buf;
	sprintf(c++, "<");
	if (N > 0) {
		sprintf(c++, "{");
		for (int i = 0; i < N; ++i) {
			sprintg(c, w[i]); while(*c)++c;
			if (i == N-1) sprintf(c++, "}");
			sprintf(c++, ",");
		}
	}
	sprintg(c, w[N]); while(*c)++c;
	sprintf(c++, ">");
	int n = c - buf;
	char *result = new char[n + 1];
	strcpy(result, buf);
	return result;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end of poly;
};

// poly type alias
typedef poly<int>    Int   ;
typedef poly<size_t> Size_t;
typedef poly<double> Double;

// Modular operation : Unit No. -> Poly Grid.
template<class T>
poly<T> operator % (T a, const poly<T>& q)
{
	int N = q.grid();
	T p[N];
	T sho = a;
	for (int i = N-1; i >= 0; --i) {
		  T qh = q[i];
		  p[i] = sho % qh;
		  sho /= qh;
	}
	return poly<T>(N, p, 0);
}

/*******************************************************************************
	This section define poly_array<T>: [N, size[N], type, data[size]]
*******************************************************************************/
template<typename T = double> class poly_array : public Int
{
template<class> friend class poly_view;
template<class> friend class poly_array_lu;

const T *data;	// data fields

void init_data(const T *dat)
{
	size_t z = this->data_size();
	T *dt = new T[z];
	memcpy(dt, dat, z * sizeof(T));
	data = dt;
}

public:
// Atom constructor
poly_array(const T *dat) : Int(Atom, 1), data(new T)
{
	T *dt = const_cast<T*>(data);
	  *dt = *dat;
}
poly_array(const T *dat, int s) : Int(Atom, s), data(new T[s])
{
	if (s == 0) throw "illegal construct, poly_array(T*,int)";
	T *dt = const_cast<T*>(data);
	memcpy(dt, dat, s * sizeof(T));
}
// Default constructor
poly_array() : data(NULL) {}
// General constructor
template <class... Args>
poly_array(const T *dat, int n, Args... cdr)   : Int(n, cdr...) { init_data(dat); }
poly_array(const T *dat, int n, int *s)        : Int(n, s)      { init_data(dat); }
poly_array(const T *dat, int n, int *s, int t) : Int(n, s, t)   { init_data(dat); }
poly_array(const T *dat, const Int& s)         : Int(s)         { init_data(dat); }
poly_array(int n, const poly_array& v)   : Int(Int(v).N-1, Int(v).w+1)
{
	size_t z = this->data_size();
	T *dt = new T[z];
	memcpy(dt, &v.data[n * z], z * sizeof(T));
	data = dt;
}
// copy constructor
poly_array(const poly_array& pa) : Int(pa) { init_data(pa.data); }
// destructor
~poly_array() { delete data; }
// substitution operator
poly_array& operator=(const poly_array& pa)
{
	if (this != &pa) {
		delete[] data;
		Int::operator=(pa);
		init_data(pa.data);
	}
	return *this;
}
poly_array& operator=(T a)
{
	size_t n = this->data_size();
	for (size_t i = 0; i < n; ++i)
		data[i] = a;
	return *this;
}
// type conversion operator
operator T* () const {return const_cast<T*>(data);}
// member function
int rows() { return Int(*this)[0]; }
int cols() { return Int(*this)[1]; }
T operator[](int pos) const
{
	return data[pos];
}
int operator()(int i) const
{
	return Int::operator[](i);
}
poly_view<T> operator[](const Int& p) const
{
	const Int& t = *this;
	int pos = t && p, N = t.N, M = p.N;
	if (pos >= t.data_size())
		throw "out of range, poly_array[Int]";
	return poly_view<T>(&data[pos], N-M, t.w + M);
}
T operator * () const { return *data; }

poly_array operator+(const poly_array& pa) const
{
	size_t n = *this; const Int& s = *this;
	size_t m = pa;
	if (n != m) throw "mismatch argument, poly_array + ";
	T result[n];
	for (size_t i = 0; i < n; ++i) result[i] = data[i] + pa.data[i];
	return poly_array(result, s);
}
poly_array operator-(const poly_array& pa) const
{
	size_t n = *this; const Int& s = *this;
	size_t m = pa;
	if (n != m) throw "mismatch argument, poly_array - ";
	T result[n];
	for (size_t i = 0; i < n; ++i) result[i] = data[i] - pa.data[i];
	return poly_array(result, s);
}
poly_array operator%(T a) const
{
	size_t n = *this; const Int& s = *this;
	T result[n];
	for (size_t i = 0; i < n; ++i) result[i] = data[i] % a;
	return poly_array(result, s);
}
char *to_s() const
{
	char buf[1024];
	char *c = buf;
	const Int& p = *this;
	char *s = p.to_s();
	sprintf(c, "{%s,", s); while(*c)++c;
	delete[] s;
	poly_to_s<T>(c, p.N, p.w, data);
	sprintf(c++, "}");
	int n = c - buf;
	c = new char[n + 1];
	strcpy(c, buf);
	return c;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end poly_array
};

template<class L, class Op, class R>
class express
{
	const L& l_; const R& r_;
public:
	express(const L& l, const R& r) : l_(l), r_(r) {}
	double operator[](int s) const { return Op::apply(l_[s], r_[s]); }
};

template<class L, class Op>
class express<L,Op,double>
{
	const L& l_; const double& r_;
public:
	express(const L& l, const double& r) : l_(l), r_(r) {}
	double operator[](int s) const { return Op::apply(l_[s], r_); }
};

template<class Op, class R>
class unary
{
	const R& r_;
public:
	unary(const R& r) : r_(r) {}
	double operator[](int s) const { return Op::apply(r_[s]); }
};

struct add { static double apply(double l, double r) { return l + r; } };
struct sub { static double apply(double l, double r) { return l - r; } };
struct mul { static double apply(double l, double r) { return l * r; } };
struct dvd { static double apply(double l, double r) { return l / r; } };
struct absolute { static double apply(double c) {return c < 0 ? -c : c;}};
struct minus { static double apply(double c) { return -c; } };

template <class R>
unary<absolute, R> abs(const R& rh)
{
	return unary<absolute, R>(rh);
}

template<typename T> void swap(poly_view<T>&, poly_view<T>&);
/*******************************************************************************
	This section define poly_view
*******************************************************************************/
template<typename T = double> class poly_view : public Int
{
template<class> friend class poly_view;
template<class> friend class poly_array_lu;

const T *data; int *cp;

public:
poly_view() : data(NULL), cp(NULL) {}
// Atom constructor
poly_view(const T *dat)        : Int(Atom, 1), data(dat), cp(NULL) {}
poly_view(const T *dat, int n) : Int(Atom, n), data(dat), cp(NULL) {}
// General constructor
poly_view(const poly_array<T>& pa)            : Int(pa),        data(pa.data), cp(NULL) {}
template <class... Args>
poly_view(const T *dat, int n, Args... cdr)   : Int(n, cdr...), data(dat), cp(NULL) {}
poly_view(const T *dat, int n, int *s)        : Int(n, s),      data(dat), cp(NULL) {}
poly_view(const T *dat, int n, int *s, int t) : Int(n, s, t),   data(dat), cp(NULL) {}
poly_view(const T *dat, const Int& s)         : Int(s),         data(dat), cp(NULL) {}
poly_view(int n, const poly_view<T>& v) : Int(Int(v).N-1, Int(v).w+1)
{
	int s = this->data_size();
	data = &v.data[n * s]; cp = NULL;
}
// copy constructor
poly_view(const poly_view& pv) : Int(pv), data(pv.data), cp(pv.cp) { cp && ++*cp; }
// destructor
~poly_view() { if (cp && --*cp == 0) { delete[] data; delete cp; }}
// substitution operator
poly_view& operator=(const poly_view& pv)
{
	if (this != &pv)
	{
		int s_size = pv.data_size();
		if (data == NULL) {
			Int::operator=(pv);
			data = new T[s_size];
			cp = new int(1);
		}
		int d_size = this->data_size();
		if ( Int::N > Atom && d_size > s_size )
		{
			int w[4] = {0}, &i = w[0], &j = w[2];
			int c = (*this)(0), n = pv(0), k(c - n), m((k + 1) / 2), o(k - m);
			Int L(&i), J(&j);
			for (i=m-1,k=c-m;i>=0;--i) { j=i+k; (*this)[J] = pv[L];      }
			for (i=n-1,k=o;  i>=0;--i) { j=i+k; (*this)[J] = pv[L];      }
			for (i=o-1,k=n;  i>=0;--i) { j=i+k; (*this)[L] = (*this)[J]; }
		}
		else if (d_size == s_size && data != pv.data) {
			T *dat = const_cast<T*>(data);
			memcpy(dat, pv.data, d_size * sizeof(T));
		}
		else if (d_size < s_size)
			throw "data size is mismatch, poly_view = ";
	}
	return *this;
}
poly_view& operator=(T a)
{
	size_t n = this->data_size();
	assert(n > 0);
	T *dat = const_cast<T*>(data);
	for (size_t i = 0; i < n; ++i) dat[i] = a;
	return *this;
}
template<class E>
poly_view& operator=(const E& e)
{
	int n = this->data_size();
	assert(n > 0);
	T *dat = const_cast<T*>(data);
	for (int i = 0; i < n; ++i) dat[i] = e[i];
	return *this;
}
// type conversion operator
operator T* () {return const_cast<T*>(data);}
// member function
int rows() {const Int& p = *this; return p[0];}
int cols() {const Int& p = *this; return p[1];}
// operators
T operator [] (int pos) const
{
	return data[pos];
}
int operator () (int i) const
{
	return Int(*this)[i];
}
poly_view operator [] (const Int& p) const
{
	const Int& t = *this;
	int pos = t && p, N = t.N, M = p.N;
	if (pos >= t.data_size())
		throw "out of range, poly_array[Int]";
	return poly_view(&data[t && p], N-M, t.w + M);
}
T operator * () const { return *data; }
template<class R> express<poly_view, add, R> operator+(const R& pv) const
{
	return express<poly_view, add, R>(*this, pv);
}
template<class R> express<poly_view, sub, R> operator-(const R& pv) const
{
	return express<poly_view, sub, R>(*this, pv);
}
template<class R> express<poly_view, mul, R> operator*(const R& pv) const
{
	return express<poly_view, mul, R>(*this, pv);
}
template<class R> express<poly_view, dvd, R> operator/(const R& pv) const
{
	return express<poly_view, dvd, R>(*this, pv);
}
unary<minus, poly_view> operator-() const
{
	return unary<minus,poly_view>(*this);
}
poly_view& operator+=(const poly_view& pv)
{
	size_t n = *this;
	for (size_t i = 0; i < n; ++i) data[i] += pv.data[i];
	return *this;
}
template<class E> poly_view& operator+=(const E& e)
{
	int n = this->data_size();
	T *dat = const_cast<T*>(data);
	for (int i = 0; i < n; ++i)
		dat[i] += e[i];
	return *this;
}
poly_view& operator-=(const poly_view& pv)
{
	size_t n = *this;
	for (size_t i = 0; i < n; ++i) data[i] -= pv.data[i];
	return *this;
}
template<class E> poly_view& operator-=(const E& e)
{
	int n = this->data_size();
	T *dat = const_cast<T*>(data);
	for (int i = 0; i < n; ++i)
		dat[i] -= e[i];
	return *this;
}
poly_view& operator*=(const poly_view& pv)
{
	size_t n = this->data_size();
	for (size_t i = 0; i < n; ++i) data[i] *= pv.data[i];
	return *this;
}
template<class E> poly_view& operator*=(const E& e)
{
	int n = this->data_size();
	for (int i = 0; i < n; ++i)
		data[i] *= e[i];
	return *this;
}
poly_view& operator/=(const poly_view& pv)
{
	size_t n = *this;
	for (size_t i = 0; i < n; ++i) data[i] /= pv.data[i];
	return *this;
}
template<class E> poly_view& operator/=(const E& e)
{
	int n = this->data_size();
	for (int i = 0; i < n; ++i)
		data[i] /= e[i];
	return *this;
}
bool operator==(const poly_view& pv)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] != pv.data[i]) break; }
	return i == n;
}
bool operator!=(const poly_view& pv)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] != pv.data[i]) break; }
	return i != n;
}
bool operator==(T a)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] != a) break; }
	return i == n;
}
bool operator!=(T a)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] != a) break; }
	return i != n;
}
bool operator<(const poly_view& pv)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] >= pv.data[i]) break; }
	return i == n;
}
bool operator>(const poly_view& pv)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] <= pv.data[i]) break; }
	return i == n;
}
bool operator<=(const poly_view& pv)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] > pv.data[i]) break; }
	return i == n;
}
bool operator>=(const poly_view& pv)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] < pv.data[i]) break; }
	return i == n;
}
bool operator<(T a)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] >= a) break; }
	return i == n;
}
bool operator>(T a)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] <= a) break; }
	return i == n;
}
bool operator<=(T a)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] > a) break; }
	return i == n;
}
bool operator>=(T a)
{
	size_t i, n = this->data_size();
	for (i = 0; i < n; ++i) { if (data[i] < a) break; }
	return i == n;
}
char *to_s() const
{
	char buf[1024];
	char *c = buf;
	const Int& p = *this;
	char *s = p.to_s();
	sprintf(c, "{%s,", s); while(*c)++c;
	sprintf(c++, "@");
	delete[] s;
	poly_to_s<T>(c, p.N, p.w, data);
	sprintf(c++, "}");
	int n = c - buf;
	c = new char[n + 1];
	strcpy(c, buf);
	return c;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end of poly_view;
};

template<typename T> void swap(poly_view<T>& p, poly_view<T>& q)
{
	size_t n = p.data_size(); T *pd = p;
	size_t m = q.data_size(); T *qd = q;
	if (n != m) throw "size mismatch, poly_view.swap";
	for (size_t i = 0; i < n; ++i) {
		T w = pd[i]; pd[i] = qd[i]; qd[i] = w;
	}
}
/*******************************************************************************
	Matrix operation
*******************************************************************************/
template <class T>
void plu_decomp(poly_view<T>& A, int *p, int &signum)
{
	if (!A.is_slate())
		throw "not Slate object in plu_decomp";
	
	int u[8]; for (int i = 0; i < 8; ++i) u[i] = 0;
	int &i = u[0], &j = u[2], &k = u[4], &pivot = u[6];
	Int L(&i), J(&j), K(&k), P(&pivot);
	int N = A.rows(), n = A.unit_size();
	T v[n], w[n];
	poly_view<T> ajj(v, n), aij(w, n);

	signum = 1; for (i=0;i<N;++i) p[i] = i;

	for (j=0;j<N-1;++j)
	{
		ajj = abs(A[J][J]);
		pivot = j;
		for (i=j+1;i<N;++i)
		{
			aij = abs(A[L][J]);
			if (aij > ajj)
			{
				ajj = aij;
				pivot = i;
			}
		}
		if (pivot != j)
		{
			auto a = A[J];
			auto b = A[P];
			swap(a, b);
			int q = p[j]; p[j] = p[pivot]; p[pivot] = q;
			signum = -signum;
		}
		if ((ajj = A[J][J]) != 0.0)
		for (i=j+1;i<N;++i)
		{
			auto au = A[L];
			auto av = A[J];
			aij = au[J] / ajj;
			au[J] = aij;
			for (k=j+1;k<N;++k) {
				au[K] -= aij * av[K];
			}
		}
		else throw "plu_decomp division by zero";
	}
}

template <class T, class S>
void plu_subst(poly_view<T>& LU, int *p, poly_view<S>& b, poly_view<S>& x)
{
	int N = LU.rows();
	int w[6] = {0}; int &i = w[0], &j = w[2], &k = w[4];
	Int L(&i), J(&j), P(&k);

	if (!(b.is_colon() || x.is_colon())) 
		throw "not Colon object in plu_subst";
	for (i=0;i<N;++i)
	{
		k = p[i]; x[L] = b[P];
	}
	for (i=1;i<N;++i)
	{
		auto sum = x[L];
		auto lu = LU[L];
		for (j=0;j<i;++j) { sum -= x[J] * lu[J]; }
		x[L] = sum;
	}
	for (i=N-1;i>=0;--i)
	{
		auto sum = x[L];
		auto lu = LU[L];
		for (j=N-1;j>i;--j) { sum -= x[J] * lu[J]; }
		x[L] = sum / lu[L];
	}
}
/*******************************************************************************
	This section define poly_array_lu
*******************************************************************************/
template <typename T> class poly_array_lu : public poly_view<T>
{

int s, *sign;	// data fields

public:
// constructor
poly_array_lu(const poly_array<T>& m) : poly_view<T>(m)
{
	size_t r = poly_view<T>::rows();
	sign = new int[r];
	plu_decomp<T>(*this, sign, s);
}
~poly_array_lu() {delete sign;}
template <typename S> poly_array<S>& solve(poly_array<S>& v)
{
	poly_array<S> *result = new poly_array<S>(v.type.N, v.type.size);
	plu_subst<T,S>(*this, sign, v, *result);
	return *result;
}
poly_array<T> invert()
{
	size_t r = poly_view<T>::rows();
	size_t c = poly_view<T>::cols();
	size_t t = poly_view<T>::unit_size();
	size_t s = r * c * t;
	T y[s], x[s];
	for (size_t i = 0; i < r; ++i)
		for (size_t j = 0; j < c; ++j)
			for (size_t k = 0; k < t; ++k)
				x[(i*c+j)*t+k] = (i == j ? 1 : 0);
	Int rct(Slate, r, c, t);
	Int ct(Colon, c, t);
	for (size_t i = 0; i < r; ++i) {
		poly_view<T> v(&y[i*c*t], ct), w(&x[i*c*t], ct);
		plu_subst<T,T>(*this, sign, w, v);
	}
	return poly_array<T>(y, rct);
}
poly_array<T> det()
{
	size_t c = poly_view<T>::cols();
	size_t t = poly_view<T>::unit_size();
	T y[t]; Int s = {0, t};
	poly_view<T> v(y, s);
	for (size_t j = 0; j < c; j++) {
		Int J(1,j,0);
		if (j == 0) v = (*this)[J][J];
		else       v *= (*this)[J][J];
	}
	return poly_array<T>(y, s);
}
template <typename S> void solve(poly_view<S>& p, poly_view<S>& a)
{
	plu_subst<T,S>(*this, sign, p, a);
}
char *to_s() const
{
	char buf[1024];
	char *c = buf;
	const poly_view<T>& v = *this;
	char *s = v.to_s();
	sprintf(c, "%s.lu", s); while(*c)++c;
	delete[] s;
	int n = c - buf;
	c = new char[n + 1];
	strcpy(c, buf);
	return c;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end poly_array_lu;
};

template <typename T>
void Kronecker(const poly_view<T>& p, const poly_view<T>& q, poly_view<T>& r)
{
	if (r.is_atom()) r = p * q;
	else
	{
		const Int& u = p;
		const Int& v = q;
		int s = u[0], t = v[0];
		for (int i = 0; i < s; i++) {
			Int L(1, i, 0);
			for (int j = 0; j < t; j++) {
				Int J(1, j, 0), K(1, i*t+j, 0);
				poly_view<T> pdr = p[L], qdr = q[J], rdr = r[K];
				Kronecker(pdr, qdr, rdr);
			}
		}
	}
}

template <typename T>
poly_array<T> operator^(const poly_array<T>& p, const poly_array<T>& q)
{
	const Int& ps = p;
	const Int& qs = q;
	Int rs = ps * qs;
	int len = rs.data_size();
	T result[len];
	poly_view<T> pv(p), qv(q), rv(result, rs);
	Kronecker(pv, qv, rv);
	return poly_array<T>(result, rs);
}

typedef poly_array<double> Poly_array;
typedef poly_view<double> Poly_view;
typedef poly_array_lu<double> Matrix_LU;

/*******************************************************************************
	comb_array<T>(T*,N,S1,...,SN,A);
	{<{S1,...,SN},A>,[(S1,...),...],...,[(SN,...),...]}
*******************************************************************************/
template<typename T = double> class comb_array : public Int
{
template <class> friend class comb_view;

const T *data;

void init_data(const T *dat)
{
	if (dat) {
		size_t s = this->comb_size();
		T *dt = new T[s];
		for (size_t i = 0; i < s; ++i) dt[i] = dat[i];
		data = dt;
	} else {
		size_t N = this->grid();
		for (size_t i = 0; i <= N; ++i) Int::w[i] = 0;
		data = NULL;
	}
}

public:
// Atom constructor
comb_array(const T *dat)        : Int(Atom, 1) { init_data(dat); }
// General constructor
template <class... Args>
comb_array(const T *dat, int n, Args... cdr)   : Int(n, cdr...) { init_data(dat); }
comb_array(const T *dat, int n, int *c)        : Int(n, c)      { init_data(dat); }
comb_array(const T *dat, int n, int *c, int t) : Int(n, c, t)   { init_data(dat); }
comb_array(const T *dat, const Int& b)         : Int(b)         { init_data(dat); }
// copy constructor
comb_array(const comb_array& ba)               : Int(ba)    { init_data(ba.data); }
// destructor
~comb_array() {delete[] data;}
// substitution operator
comb_array& operator=(const comb_array& ba)
{
	if (this != &ba) {
		delete data;
		Int::operator=(ba);
		init_data(ba.data);
	}
	return *this;
}
const T* operator[](int i) const
{
	int s = this->cast_size(i);
	return &data[s];
}
int operator()(int i) const
{
	return Int::operator[](i);
}
char *to_s() const
{
	char buf[1024];
	char *c = buf;
	const Int& p = *this;
	char *s = p.to_s();
	sprintf(c, "{%s,", s); while(*c)++c;
	delete[] s;
	int N = p.grid(), n = 0, A = p.unit_size();
	if (N > 0) for (int i = 0; i < N; ++i) {
		sprintf(c++, "[");
		int l = p[i];
		for (int j = 0; j < l; ++j) {
			const T *pos = &data[n+j*A];
			if (A > 1) sprintf(c++, "(");
			for (int a = 0; a < A; ++a) {
				sprintg(c, pos[a]); while(*c)++c;
				if (a < A-1) {sprintf(c++, ",");}
			}
			if (A > 1) {sprintf(c++, ")");}
			if (j < l-1) {sprintf(c++, ",");}
		}
		sprintf(c++, "]");
		if (i < N-1) {sprintf(c++, ",");}
		n += l * A;
	} else
		poly_to_s<T>(c, N, p.w, data);
	sprintf(c++, "}");
	n = c - buf;
	c = new char[n + 1];
	strcpy(c, buf);
	return c;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end of comb_array;
};

template<typename T = double> class comb_view : public Int
{

const T *data;

public:
// Atom constructor
comb_view(const T *dat)        : Int(Atom, 1), data(dat) {}
comb_view(const T *dat, int n) : Int(Atom, n), data(dat) {}
// General constructor
template <class... Args>
comb_view(const T *dat, int n, Args... cdr)   : Int(n, cdr...), data(dat) {}
comb_view(const T *dat, int n, int *c)        : Int(n, c),      data(dat) {}
comb_view(const T *dat, int n, int *c, int t) : Int(n, c, t),   data(dat) {}
comb_view(const T *dat, const Int& b)         : Int(b),         data(dat) {}
comb_view(const comb_array<T>& ca)            : Int(ca),    data(ca.data) {}
// copy constructor
comb_view(const comb_view& ba) : Int(ba), data(ba.data) {}
// destructor
// substitution operator
comb_view& operator=(const comb_view& ba)
{
	if (this != &ba) {
		int s = this->comb_size();
		int t = ba.comb_size();
		if (s != t) throw "size mismatch, comb_view = ";
		for (int i = 0; i < s; ++i) data[i] = ba.data[i];
	}
	return *this;
}
const T* operator[](int i) const
{
	int n = this->is_atom() ? *Int::w : Int::N;
	if (i < 0 || i >= n) throw "out of range, comb_view[]";
	int s = this->cast_size(i);
	return &data[s];
}
int operator()(int i) const
{
	if (i < 0 || i > Int::N) throw "out of range, comb_view()";
	return Int::operator[](i);
}
char *to_s() const
{
	char buf[1024];
	char *c = buf;
	const Int& p = *this;
	char *s = p.to_s();
	sprintf(c, "{%s,", s); while(*c)++c;
	delete[] s;
	int N = p.grid(), n = 0, A = p.unit_size();
	if (N > 0) for (int i = 0; i < N; ++i) {
		sprintf(c, "@["); while(*c)++c;
		int l = p[i];
		for (int j = 0; j < l; ++j) {
			const T *pos = &data[n+j*A];
			if (A > 1) {sprintf(c++, "(");}
			for (int a = 0; a < A; ++a) {
				sprintg(c, pos[a]); while(*c)++c;
				if (a < A-1) {sprintf(c++, ",");}
			}
			if (A > 1) {sprintf(c++, ")");}
			if (j < l-1) {sprintf(c++, ",");}
		}
		sprintf(c++, "]");
		if (i < N-1) {sprintf(c++, ",");}
		n += l * A;
	} else {
		sprintf(c++, "@");
		poly_to_s<T>(c, N, p.w, data);
	}
	sprintf(c++, "}");
	n = c - buf;
	c = new char[n + 1];
	strcpy(c, buf);
	return c;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end of comb_view;
};

template<typename T> class comp_array : public Int
{
	T *data;

public:
// mono constructor
comp_array(T *dat, int o, int s, int k = 1) : Int(2, o, s, k), data(new T[s*k])
{
	size_t sk = s * k;
	for (size_t i = 0; i < sk; ++i) data[i] = dat[i];
}
// General constructor
comp_array(T *dat, int n, int *s) : Int(n, s)
{
	int N = this->grid();
	size_t sk = this->cast_size(N, N/2);
	data = new T[sk];
	memcpy(data, dat, sk * sizeof(T));
}
comp_array(T *dat, int n, int *s, int k) : Int(n, s, k)
{
	int N = this->grid();
	size_t sk = this->cast_size(N, N/2);
	data = new T[sk];
	memcpy(data, dat, sk * sizeof(T));
}
// copy constructor
comp_array(const comp_array& pa) : Int(pa)
{
	int N = this->grid();
	size_t sk = this->cast_size(N, N/2);
	data = new T[sk];
	memcpy(data, pa.data, sk * sizeof(T));
}
// destructor
~comp_array() {delete data;}
// substitution operator
comp_array& operator=(const comp_array& pa)
{
	if (this != &pa) {
		delete data;
		Int::operator=(pa);
		int N = this->grid();
		size_t sk = this->cast_size(N, N/2);
		data = new T[sk];
		memcpy(data, pa.data, sk * sizeof(T));
	}
	return *this;
}
// member function
Int offset() const
{
	const Int& b = *this;
	int N = b.grid();
	return Int(N/2, b.w, 0);
}
Int size() const
{
	const Int& b = *this;
	int N = b.grid();
	return Int(N/2, &b.w[N/2], b.unit_size());
}
comp_view<T> operator[](const Int& p) const
{
	const Int& b = *this;
	int N = b.grid(), K = b.unit_size();
	T **nk = new T*[N/2];
	for (int n = 0, i = 0; i < N/2; ++i) {
		int s = b[i+N/2];
		int k = n + p[i];
		nk[i] = &data[k];
		n += s * K;
	}
	return comp_view<T>(nk, N/2, &(b.w[N/2]), p.w[N/2]);
}
T* operator[](int i) const
{
	const Int& b = *this;
	int N = b.grid();
	return &data[b.cast_size(N/2 + i, N/2)];
}
int operator() (int i) const
{
	return Int::w[i];
}
operator T* () {return data;}
comp_array operator & (const comp_array& pa)
{
	const Int& bt = *this;
	const Int& pt = pa;
	int Nt = bt.grid(), Np = pt.grid();
	int s[Nt + Np + 1], *sp = s;
	for (int i = 0; i < Nt/2; ++i) *sp++ = bt[i];
	for (int i = 0; i < Np/2; ++i) *sp++ = pt[i];
	for (int i =Nt/2; i < Nt; ++i) *sp++ = bt[i];
	for (int i =Np/2; i < Np; ++i) *sp++ = pt[i];
	size_t skt = bt.cast_size(Nt, Nt/2);
	size_t skp = pt.cast_size(Np, Np/2);
	T dt[skt + skp], *dp = dt;
	for (size_t i = 0; i < skt; ++i) *dp++ = data[i];
	for (size_t i = 0; i < skp; ++i) *dp++ = pa.data[i];
	return comp_array(dt, Nt + Np, s, bt.unit_size());
}
char *to_s() const
{
	char buf[1024];
	char *c = buf;
	const Int& p = *this;
	int N = p.grid(), K = p.unit_size();
	sprintf(c, "{<"); while(*c)++c;
	if (N > 0) {
		sprintf(c++, "{");
		for (int i = 0; i < N/2; ++i) {
			sprintf(c, "%d", p[i]); while(*c)++c;
			if (i == N/2-1) sprintf(c++, "}");
			sprintf(c++, ",");
		}
		sprintf(c++, "{");
		for (int i = N/2; i < N; ++i) {
			sprintf(c, "%d", p[i]); while(*c)++c;
			if (i == N-1) sprintf(c++, "}");
			sprintf(c++, ",");
		}
	}
	sprintf(c, "%d>,",p[N]); while(*c)++c;
	int n = 0;
	for (int i = 0; i < N/2; ++i) {
		int z = p[i + N/2];
		sprintf(c++, "[");
		for (int k = 0; k < K; ++k) {
			sprintf(c++, "(");
			for (int j = 0; j < z; ++j) {
				sprintg(c, data[n+z*k+j]); while(*c)++c;
				if (j < z-1) {sprintf(c++, ",");}
			}
			sprintf(c++, ")");
			if (k < K-1) {sprintf(c++, ",");}
		}
		sprintf(c++, "]");
		if (i < N/2-1) {sprintf(c++, ",");}
		n += z * K;
	}
	sprintf(c++, "}");
	n = c - buf;
	c = new char[n + 1];
	strcpy(c, buf);
	return c;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end of comp_array;
};

template<typename T> class comp_view : public Int
{
	T **data;
public:
template <class... Args>
comp_view(T **dat, int n, int car, Args... cdr) : Int(n, car, cdr...), data(dat) {}
comp_view(T **dat, int n, const int *c        ) : Int(n, c          ), data(dat) {}
comp_view(T **dat, int n, const int *c, int t ) : Int(n, c, t       ), data(dat) {}
comp_view(T **dat, const Int& b               ) : Int(b             ), data(dat) {}
~comp_view() {delete[] data;}
T operator * () const
{
	return coeff(data, Int::w, this->grid());
}
T operator[](int j) const
{
	return coeff(data, Int::w, this->grid(), j);
}
T operator[](T *ds) const
{
	return total_derive(data, Int::w, this->grid(), this->unit_size()-1, ds);
}
T operator () (T dh, int j = 0) const
{
	const Int& p = *this;
	int k = p.unit_size();
	T sum = 0.0;
	T *ds = data[j];
	int c = p[j] * k;
	for (int i = k; i > 0; --i)
		sum = (sum + ds[c -= p[j]]) * dh / (T)i;
	return sum;
}

T taylor(T dh, int j = 0) const
{
	const Int& p = *this;
	int k = p.unit_size();
	T sum = 0.0;
	T *ds = data[j];
	int c = p[j] * k;
	for (int i = k; i > 0; --i)
		sum = sum * dh / (T)i + ds[c -= p[j]];
	return sum;
}
char *to_s() const
{
	char buf[1024];
	char *c = buf;
	const Int& p = *this;
	int N = p.grid();
	int K = p.unit_size();
	char *s = p.to_s();
	sprintf(c, "{%s", s); while(*c)++c;
	delete[] s;
	sprintf(c, "@["); while(*c)++c;
	for (int i = 0; i < N; ++i) {
		if (K > 1) sprintf(c++, "[");
		for (int j = 0; j < K; ++j) {
			sprintg(c, data[i][j*p[i]]); while(*c)++c;
			if (j < K-1) sprintf(c++, ",");
		}
		if (K > 1) sprintf(c++, "]");
		if ( i < N-1) sprintf(c++, ",");
	}
	sprintf(c++, "]");
	int n = c - buf;
	c = new char[n + 1];
	strcpy(c, buf);
	return c;
}
void print(const char *end = NULL) const {
	char *s = this->to_s();
	char f[32] = {"%s"};
	char *c = f; while(*c)++c;
	if (end) sprintf(c, "%s", end);
	fprintf(stdout, f, s);
	delete[] s;
}
// end of comp_view;
};

template <typename T>
bool xvalue(int &i, T &x, const T *X, T &t, int &l, int L, int &d, int p, bool F)
{
	bool f = false;
	if (F) if (d++ == p) {
		if (++l == L)  {
			f = true;
			l = 1;
		}	t = (X ? (X[l] - X[l-1]) : 1.0)/p;
		d = (l == 1) ? 0 : 1;
	}
	i = d + (l - 1) * p;
	x = (X ? X[l - 1] : l - 1) + d * t;
	return f;
}

template <typename T> class poly_index
{
	int N, *index, *Le, *Lc, *Nc, Dx;
	T *value, *Dt; const T **X; Int *e;
public:
poly_index(const comb_array<T>& x, int Dp)
{
	const Int& end = x; N = end.grid(); int Np[N]; 
	index = new int[N*4]; Le = index + N; Lc = Le + N; Nc = Lc + N;
	value = new T[N*2]; Dt = value + N; X = new const T*[N];
	for (int i = 0; i < N; i++) {
		Np[i] =(end[i] - 1) * Dp + 1;
		Le[i] = end[i]; X[i] = x[i];
		Lc[i] = 1; Nc[i] = 0;
		Dt[i] = (X[i] ? (X[i][1] - X[i][0]) : 1.0) / Dp;
		xvalue(index[i], value[i], X[i], Dt[i], Lc[i], Le[i], Nc[i], Dp, false); 
	}	Dx = Dp;
	e = new Int(N, Np, 1);
}
~poly_index() {delete e; delete X; delete value; delete index;}
poly_index& operator++()
{
	bool f = true;
	for (int i = N-1; i >= 0; --i)
		f = xvalue(index[i], value[i], X[i], Dt[i], Lc[i], Le[i], Nc[i], Dx, f); 
	return *this;
}
int index_size()	// return Npx;
{
	return e->grid_size();
}
int index_value()	// convert index to number;
{
	Int i(N, index, 0); 
	return (*e) && i;
}
Int get_index() {return Int(N, index, 0);}
poly<T> get_value() {return poly<T>(N, value, 0);}
// end of poly_index
};

#endif
