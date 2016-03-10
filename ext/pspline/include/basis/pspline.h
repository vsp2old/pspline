#ifndef _PSPLINE_5_H_INCLUDED_
#define _PSPLINE_5_H_INCLUDED_

template <typename T>
void poly_new_subst(int n, T **w)
{
	for (int i = 0; i < n; ++i) w[i] = NULL;
}

template <typename T, class... Args>
void poly_new_subst(int n, T **w, const T *car, const Args... cdr)
{
	w[0] = new T(*car); if (n > 1) poly_new_subst(n-1, &w[1], cdr...);
}

/*******************************************************************************
	define poly_spline[base_spline1*,...base_splineN*]
*******************************************************************************/
template <typename S = double> class poly_spline
{

int N; base_spline<S> **base;

public:
struct values_ {
	Int imax, icox, rank, jisu, shuki;
	carray_view<S> *knots;
	values_(poly_spline<S>& own, int as) :
		imax(own.Imax(as)), icox(own.Icox(as)),
		rank(own.Rank()  ), jisu(own.Jisu()  ), shuki(own.Shuki())
		{
			size_t n = own.Grid();
			knots = new carray_view<S>(n, own.Maxq(), NULL);
			for (size_t i = 0; i < n; i++)
				(*knots)[i] = own[i].Knots();
		}
	~values_() { delete knots; }
} *values;

// constructor
poly_spline() : base(NULL), values(NULL) {}
poly_spline(const comb_array<S>& x, const Int& n, const Int& j, const Int& s) : N(x.grid())
{
	base = new base_spline<S>*[N];
	for (int i = 0; i < N; i++)
		base[i] = new base_spline<S>(n[i], x[i], j[i], -s[i]%2);
	values = new values_(*this, n.unit_size());
}
template <class... Args>
poly_spline(int k, const base_spline<S> *bsc, Args... bsd)
 : N(sizeof...(bsd)+1), base(new base_spline<S>*[sizeof...(bsd)+1])
{
	base[0] = new base_spline<S>(*bsc);
	if (N > 1) poly_new_subst(N-1, &base[1], bsd...);
	values = new values_(*this, k);
}
poly_spline(int k, const poly_spline<S> *ps, const base_spline<S> *bs)
 : N(ps->Grid()+1), base(new base_spline<S>*[ps->Grid()+1])
{
	for (int i = 0; i < N-1; i++) {
		const base_spline<S>& b = (*ps)[i];
		base[i] = new base_spline<S>(b);
	}
	base[N-1] = new base_spline<S>(*bs);
	values = new values_(*this, k);
}
poly_spline(int k, const base_spline<S> *bs, const poly_spline<S> *ps)
 : N(ps->Grid()+1), base(new base_spline<S>*[ps->Grid()+1])
{
	base[0] = new base_spline<S>(*bs);
	for (int i = 1; i < N; i++) {
		const base_spline<S>& b = (*ps)[i-1];
		base[i] = new base_spline<S>(b);
	}
	values = new values_(*this, k);
}
poly_spline(int k, const base_spline<S>& bsc)
 : N(1), base(new base_spline<S>*[1])
{
	base[0] = new base_spline<S>(bsc);
	values = new values_(*this, k);
}
// copy constructor
poly_spline(const poly_spline& p) : N(p.N), base(new base_spline<S>*[p.N])
{
	for (int i = 0; i < N; i++)
		base[i] = new base_spline<S>(*p.base[i]);
	values = new values_(*this, p.values->imax.unit_size());
}
// destructor
~poly_spline() {delete[] base; delete values;}
// substitution operator
poly_spline& operator=(const poly_spline& p)
{
	if (this != &p) {
		delete[] base;
		N = p.N; base = new base_spline<S>*[p.N];
		for (int i = 0; i < N; i++) {
			base[i] = new base_spline<S>(*p.base[i]);
		}	delete values;
		values = new values_(*this, p.values->imax.unit_size());
	}
	return *this;
}

marray<S> gyoretsu(base_spline<S> *base, int a, const S *x)
{
	int imax = base->Imax();
	marray<S> bd(imax, imax);
	base->gyoretu(bd, a, x);
	return bd;
}

marray<S> gyoretu(const Int& a, const comb_array<S>& x)
{
	marray<S> res = gyoretsu(base[0], a[0], x[0]);
	for (int i = 1; i < N; ++i) {
		res = res ^ gyoretsu(base[i], a[i], x[i]);
	}
	return res;
}

comp_array<S> basis(const poly<S>& t, const int *jb = NULL) const
{
	int s[N*2];
	int size = this->values->rank.comb_size();
	S data[size], *p = data;
	for (int i = 0; i < N; i++) {
		parray<S> bc = base[i]->basis(t[i], jb ? jb[i] : 0);
		s[i] = bc(0); s[i+N] = bc(1);
		S *dat = *bc;
		for (int j = 0; j < bc(1); ++j) *p++ = dat[j];
	}
	return comp_array<S>(data, N*2, s, 1);
}

comp_array<S> bases(const poly<S>& t, int jbn = 0) const
{
	int s[N*2], K = jbn + 1;
	int size = this->values->rank.comb_size();
	S data[size * K], *p = data;
	for (int i = 0; i < N; i++) {
		parray<S> bc = base[i]->bases(t[i], jbn);
		s[i] = bc(0); s[i+N] = bc(1);
		int sz = bc(1) * K;
		S *dat = *bc;
		for (int j = 0; j < sz; ++j) *p++ = dat[j];
	}
	return comp_array<S>(data, N*2, s, K);
}

template <typename T>
void sekibun(poly_view<T>& r, const poly_view<T>& clp, const poly<int>& h, const poly<S>& x) const
{
	int offset[N+1], size[N+1];
	Size_t z = this->values->rank;
	S **sp = create_carray<S>(N, z);
	for (int i = 0; i < N; ++i) {
		S ql = (*base[i])[h[i]];
		S qm = (*base[i])[h[i] + 1];
		S qh = (x[i] < qm ? x[i] : qm) - ql;
		parray<S> pc = base[i]->bases(ql, base[i]->Jisu());
		offset[i] = pc(0); size[i] = pc(1);
		for (int j = 0; j < size[i]; j++) sp[i][j] = pc.saylor(j, qh);
	}	offset[N] = 0; size[N] = 0;
	Int o = Int(offset,N);
	Int s = Int(size,  N);
	int end = s.grid_size();
	for (int i = 0; i < end; ++i) {
		Int k = i % s, j = k + o;
		r += clp[j] * [](int N, T **p, Int& k) -> T {
			T u = 1.0;
			for (int i = 0; i < N; ++i) u *= p[i][k[i]];
			return u;
		} (N, sp, k);
	}
	FREE(sp);
}

base_spline<S>& operator[](int i) const {return *base[i];}

int Grid() const {return N;}

Int Kset(const poly<S>& x) const
{
	int v[N];
	for (int i = 0; i < N; i++) v[i] = base[i]->Kset(x[i]);
	return Int(N, v, 0);
}

private:

Int Imax(int as = 1) const
{
	int v[N];
	for (int i = 0; i < N; i++) v[i] = base[i]->Imax();
	return Int(N, v, as);
}
Int Icox(int as = 1) const
{
	int v[N];
	for (int i = 0; i < N; i++) v[i] = base[i]->Icox();
	return Int(N, v, as);
}
Int Rank() const
{
	int v[N];
	for (int i = 0; i < N; i++) v[i] = base[i]->Rank();
	return Int(N, v, 1);
}
Int Jisu() const
{
	int v[N];
	for (int i = 0; i < N; i++) v[i] = base[i]->Jisu();
	return Int(N, v, 0);
}
Int Shuki() const
{
	int v[N];
	for (int i = 0; i < N; i++) v[i] = base[i]->Shuki();
	return Int(N, v, 0);
}
size_t *Maxq() const
{
	size_t *vp = new size_t[N];
	for (int i = 0; i < N; i++) vp[i] = base[i]->Maxq();
	return vp;
}
// end poly_spline
};

template <class T, class S> class bspline;
/*******************************************************************************
	Ｎ変数Ｋ次元パラメトリック ＆ リーゼンフェルト スプライン補間
	     int N : number of parameters (Grid)
	typename S : poly<S,N> = (s, t, u, ...)
	     int K : number of axis (Atom)
	typename T : poly<T,K> = (x, y, z, ...)
--------------------------------------------------------------------------------
	パラメトリックスプライン      （曲面補間）  s = {0,...}
	パラメトリック周期スプライン  （閉曲面補間）s = {1,...}
	リーゼンフェルトスプライン    （曲面補間）  s = {2,...}
	リーゼンフェルト周期スプライン（閉曲面補間）s = {3,...}
*******************************************************************************/
template <typename T, typename S = double>
class pspline : public poly_spline<S>
{

	int K; T *clp;

  public:
	// constructor
	pspline():clp(NULL){}
	pspline(const poly_array<T>&, const comb_array<S>&, const Int&, const Int&);
	template <class... Args>
	pspline(int k, T *v, Args... bset)
	 : poly_spline<S>(k, bset...), K(k), clp(v) {}
	pspline(const bspline<T,S>& bs)
	 : poly_spline<S>(bs.Unit_size(), base_spline<S>(bs)), K(bs.Unit_size())
	{
		const Int& c = poly_spline<S>::values->icox;
		int s = c.data_size();
		clp = new T[s];
		memcpy(clp, (T*)bs, s * sizeof(T));
	}
	pspline(const pspline&);
	// destructor
	~pspline() { delete[] clp; }
	// operator
	pspline& operator=(const pspline&);
	poly_array<T> operator[](S t) const { return (*this)(t); }
	poly_array<T> operator()(S, int = 0) const;
	poly_array<T> operator[](const poly<S>& p) const {return (*this)(p);}
	poly_array<T> operator()(const poly<S>&, const int* = NULL) const;
	poly_array<T> operator()(const poly<S>&, int, T*) const;
	poly_array<T> sekibun(const poly<S>&) const;
	int Unit_size() const {return K;}
	operator T* () const {return clp;}
};

template <typename T, typename S>
pspline<T,S>::pspline(const poly_array<T>& p, const comb_array<S>& x, const Int& j, const Int& s)
 : poly_spline<S>(x, Int(p), j, s), K(p.unit_size())
{
	const Int& a = poly_spline<S>::values->imax;
	const Int& c = poly_spline<S>::values->icox;
	T *alp = clp = new T[c.data_size()];
	poly_view<T> alv(alp, a), clv(clp, c);
	if (s < 2) {
		marray<S> bd = poly_spline<S>::gyoretu(a, x);
		marray_lu_view<S> lu(bd);
		lu.solve((T*)p, alp, K);
	} else
		memcpy(alp, (T*)p, p.data_size() * sizeof(T));
	clv = alv;
}
// copy constructor
template <typename T, typename S>
pspline<T,S>::pspline(const pspline& ps) : poly_spline<S>(ps), K(ps.K)
{
	const Int& c = poly_spline<S>::values->icox;
	int len = c.data_size();
	clp = new T[len];
	memcpy(clp, (T*)ps, len * sizeof(T));
}
// substitution operator
template <typename T, typename S>
pspline<T,S>& pspline<T,S>::operator=(const pspline& ps)
{
	if (this != &ps) {
		delete clp;
		poly_spline<S>::operator=(ps); K = ps.K;
		const Int& c = poly_spline<S>::values->icox;
		int len = c.data_size();
		clp = new T[len];
		memcpy(clp, (T*)ps, len * sizeof(T));
	}
	return *this;
}

template <typename T, typename S>
poly_array<T> pspline<T,S>::operator()(const poly<S>& t, const int *jb) const
{
	if (jb) {
		Int jbn(t.grid(), jb, 0);
		if (!(poly_spline<S>::values->jisu >= jbn))
			throw "out of range, pspline()";
	}
	T result[K]; for (int i = 0; i < K; ++i) result[i] = 0;
	poly_view<T> res(result, K);
	comp_array<S> pc = poly_spline<S>::basis(t, jb);
	const Int& c = poly_spline<S>::values->icox;
	poly_view<T> clv(clp, c);
	Int offset = pc.offset(), size = pc.size();
	int s = size.grid_size();
	for (int i = 0; i < s; ++i)
	{
		Int k = i % size, l = k + offset;
		res += clv[l] * pc[k = 1][0];
	}
	return poly_array<T>(result, K);
}

template <typename T, typename S>
poly_array<T> pspline<T,S>::operator()(const poly<S>& t, int jbn, T *ds) const
{
	if (!(poly_spline<S>::values->jisu >= jbn))
		throw "out of range, pspline()";
	T result[K]; for (int i = 0; i < K; ++i) result[i] = 0;
	poly_view<T> res(result, K);
	comp_array<S> pc = poly_spline<S>::bases(t, jbn);
	const Int& c = poly_spline<S>::values->icox;
	poly_view<T> clv(clp, c);
	Int offset = pc.offset(), size = pc.size();
	int s = size.grid_size();
	for (int i = 0; i < s; ++i)
	{
		Int k = i % size, l = k + offset;
		res += clv[l] * pc[k = jbn + 1][ds];
	}
	return poly_array<T>(result, K);
}

template <typename T, typename S>
poly_array<T> pspline<T,S>::sekibun(const poly<S>& t) const
{
	const Int& c = poly_spline<S>::values->icox;
	const Int& j = poly_spline<S>::values->jisu;
	Int n = poly_spline<S>::Kset(t) - j + 1;
	int Imax = n.grid_size();
	poly_view<T> clw(clp, c);
	int a = c.unit_size();
	T result[a]; for (int i = 0; i < a; i++) result[i] = 0;
	poly_view<T> R(result, a);
	for (int i = 0; i < Imax; ++i) {
		Int h = i % n + j;
		poly_spline<S>::sekibun(R, clw, h, t);
	}
	return poly_array<T>(result, a);
}

template <typename T, typename S>
poly_array<T> pspline<T,S>::operator()(S t, int b) const
{
	int N = this->Grid();
	S w[N+1]; for (int i = 0; i < N; i++) w[i] = t; w[N] = 0;
	poly<S> u(w, N);
	if (!(poly_spline<S>::values->jisu >= b))
		throw "out of range, pspline()";
	T result[K]; for (int i = 0; i < K; ++i) result[i] = 0;
	poly_view<T> res(result, K);
	comp_array<S> pc = poly_spline<S>::bases(u, b);
	const Int& c = poly_spline<S>::values->icox;
	poly_view<T> clv(clp, c);
	Int offset = pc.offset(), size = pc.size();
	int s = size.grid_size();
	for (int i = 0; i < s; ++i)
	{
		Int k = i % size, l = k + offset;
		res += clv[l] * pc[k = b + 1][b];
	}
	return poly_array<T>(result, K);
}

//	１変数Ｎ次元パラメトリックスプライン
typedef pspline<double,double> Pspline;
//	１変数Ｎ次元パラメトリックスプライン
typedef pspline<double,double> Pspline2;

/*******************************************************************************
	１変数多次元パラメトリックスプライン
*******************************************************************************/
template <class T = double, class S = T> class bspline : public base_spline<S>
{
	int K; T *clp;

  public:
	// constructor
	bspline(const poly_array<T>&, int, S*, int, int = 0, int* = NULL);
	bspline(const bspline&);
	// destructor
	~bspline() { delete clp; }
	// operator
	bspline& operator=(const bspline&);
	poly_array<T> operator[](S t) const {return (*this)(t);}
	poly_array<T> operator()(S, int = 0) const;
	poly_array<T> sekibun(S) const;
	T line_element(S) const;
	bspline *line_integral(T(*)(poly_array<S>&), int, S*, S = 0, T = 0) const;
	int Unit_size() const {return K;}
	operator T* () const {return clp;}
};
/*******************************************************************************
	１変数パラメトリックスプライン コンストラクタ
--------------------------------------------------------------------------------
	パラメトリックスプライン      （曲線補間）  s = 0
	パラメトリック周期スプライン  （閉曲線補間）s = 1
	リーゼンフェルトスプライン    （曲線補間）  s = 2
	リーゼンフェルト周期スプライン（閉曲線補間）s = 3
*******************************************************************************/
template <class T, class S>
bspline<T,S>::bspline(
	const poly_array<T>& p, // Vector Array<1>, [P0,...,Pn-1];
			  int n, S*  x, // Variables,       [X0,...,Xn-1];
					int  j, // Dimension
					int  s, // Type parameter
					int *d	// Differential rank, [D0,...,Dn-1];
)
: base_spline<S>(n, x, j, p(0) - n - (s%2)), K(p.unit_size())
{
	int a = base_spline<S>::imax;
	int c = base_spline<S>::icox;
	T *alp = clp = new T[c * K];
	if (s < 2) {
		marray<S> bd(a, a);
		base_spline<S>::gyoretu((S**)bd, a, x, d);	//係数行列を作成
		marray_lu_view<S> lu(bd);
		lu.solve((T*)p, alp, K);
	} else
		memcpy(alp, (T*)p, a * K * sizeof(T));
	poly_view<T> alv(alp, Colon, a, K);
	poly_view<T> clv(clp, Colon, c, K);
	clv = alv;
}

template <class T, class S>
bspline<T,S>::bspline(const bspline& s) : base_spline<S>(s), K(s.K)
{
	int c = base_spline<S>::icox;
	clp = new T[c * K];
	memcpy(clp, (T*)s, c * K * sizeof(T));
}

template <class T, class S>
bspline<T,S>& bspline<T,S>::operator=(const bspline& s)
{
	if (this != &s) {
		base_spline<S>::operator=(s);
		K = s.K;
		delete clp;
		int c = base_spline<S>::icox;
		clp = new T[c * K];
		memcpy(clp, (T*)s, c * K * sizeof(T));
	}
	return *this;
}

template <typename T, typename S>
void sekibun(const bspline<T,S>& p, S t, const poly_view<T>& alp, int Jsk, poly_view<T>& sek)
{
	const base_spline<S>& base = p;
	const Int& size = alp;
	int ks, kset, k = base.Rank() + Jsk;
	S b[k], *q = base.Knots();
	ks = deboor<S>(base.Rank(), q, base.Icox(), t, b, kset, k);
	poly_view<T> *a[2], *w = new poly_view<T>[2*(kset+1)];
	a[0] = w; a[1] = &w[kset+1];
	for (int i = 0; i <= kset; ++i) {
		Int k = i % size;
		a[Jsk%2][i] = alp[k]; a[(Jsk+1)%2][i] = sek;
	}
	while (Jsk) {
		poly_view<T> *a0 = a[Jsk%2], *a1 = a[(Jsk+1)%2];
		for (int i = 0; i <= kset; ++i) {
			if (i == 0) a1[i] = 0.0;
			else		a1[i] = a1[i-1];
			a1[i] += a0[i] * ((q[i+(k-Jsk)] - q[i]) / S(k-Jsk));
		}
		Jsk--;
	}
	for (int i = ks > 0 ? ks : 0; i <= kset; ++i)
		sek += a[0][i] * b[i-ks];
	delete[] w;
}

//
//	１変数パラメトリックスプライン関数値の計算（曲線補間）
//

template <class T, class S>
poly_array<T> bspline<T,S>::operator()(S t, int b) const
{
	if (b > base_spline<S>::jisu)
		throw "out of range, bspline value";
	T result[K]; for (int i = 0; i < K; ++i) result[i] = 0;
	poly_view<T> res(result, size_t(K));
	Int c(Colon, base_spline<S>::icox, K);
	poly_view<T> clv(clp, c);
	if (b >= 0) {
		parray<S> bc = base_spline<S>::basis(t, b);
		Int offset(1, bc(0), 0), size(1, bc(1), 0);
		int s = size.grid_size();
		for (int i = 0; i < s; ++i)
		{
			Int k = i % size, l = k + offset;
			res += clv[l] * bc[0][i];
		}
	} else ::sekibun<T,S>(*this, t, clv, -b, res);
	return poly_array<T>(result, K);
}

template <class T, class S>
poly_array<T> bspline<T,S>::sekibun(S t) const
{
	T result[K]; for (int i = 0; i < K; ++i) result[i] = 0.0;
	base_spline<S>::sekibun(t, clp, K, result);
	return poly_array<T>(result, K);
}

template <class T, class S>
T bspline<T,S>::line_element(S t) const
{
	poly_array<T> r = (*this)(t, 1);
	T s = 0;
	for (int i = 0; i < K; ++i) s += r[i] * r[i];
	return sqrt(s);
}

template <class T, class S>
bspline<T,S> *bspline<T,S>::line_integral(T(*F)(poly_array<S>&), int n, S *u, S x, T y) const
{
	int N = n+1, D[N]; S X[N]; T Y[N];
	for (int i = 0; i < n; i++) D[i] = 1; D[n] = 0;
	for (int i = 0; i < n; i++) {
		poly_array<S> U = (*this)[u[i]];
		X[i] = u[i];
		Y[i] = F(U) * this->line_element(u[i]);
	}	X[n] = x;
		Y[n] = y;
	poly_array<T> Ya(Y, 1, N, 1);
	return new bspline<T,S>(Ya, n, X, base_spline<S>::jisu, 0, D);
}

template <class T, class S>
pspline<T,S> operator * (const bspline<T,S>& bs1, const bspline<T,S>& bs2)
{
	int k1 = bs1.Unit_size(), k2 = bs2.Unit_size();
	assert(k1 == k2);
	int c1 = bs1.Icox(), c2 = bs2.Icox();
	T *cp1 = bs1, *cp2 = bs2, *cp3 = new T[c1*c2*k1];
	for (int i = 0; i < k1; i++) {
		varray_view<T> v1(cp1+i, c1, k1), v2(cp2+i, c2, k2), v3(cp3+i, c1*c2, k1);
		kronecker(v1, v2, v3);
	}
	const base_spline<S>& bt1 = bs1;
	const base_spline<S>& bt2 = bs2;
	return pspline<T,S>(k1, cp3, &bt1, &bt2);
}

template <class T, class S>
pspline<T,S> operator * (const pspline<T,S>& ps, const bspline<T,S>& bs2)
{
	int k2 = bs2.Unit_size(), c2 = bs2.Icox();
	const poly_spline<S>& ps1 = ps;
	const Int& c = ps1.values->icox;
	int k1 = c.unit_size(), c1 = c.grid_size();
	assert(k1 == k2);
	T *cp1 = ps, *cp2 = bs2, *cp3 = new T[c1*c2*k2];
	for (int i = 0; i < k1; i++) {
		varray_view<T> v1(cp1+i, c1, k1), v2(cp2+i, c2, k2), v3(cp3+i, c1*c2, k1);
		kronecker(v1, v2, v3);
	}
	const base_spline<S>& bt2 = bs2;
	return pspline<T,S>(k2, cp3, &ps1, &bt2);
}

template <class T, class S>
pspline<T,S> operator * (const bspline<T,S>& bs1, const pspline<T,S>& ps)
{
	int k1 = bs1.Unit_size(), c1 = bs1.Icox();
	const poly_spline<S>& ps2 = ps;
	const Int& c = ps2.values->icox;
	int k2 = c.unit_size(), c2 = c.grid_size();
	assert(k1 == k2);
	T *cp1 = bs1, *cp2 = ps, *cp3 = new T[c1*c2*k1];
	for (int i = 0; i < k1; i++) {
		varray_view<T> v1(cp1+i, c1, k1), v2(cp2+i, c2, k2), v3(cp3+i, c1*c2, k1);
		kronecker(v1, v2, v3);
	}
	const base_spline<S>& bt1 = bs1;
	return pspline<T,S>(k1, cp3, &bt1, &ps2);
}
//
//	misc
//
// 偏微分プロット
template <class T, class S = double>
int plot(const pspline<T,S>&, const comb_array<S>&, S**, T**, int, const int* = NULL);
// 全微分プロット
template <class T, class S = double>
int plot(const pspline<T,S>&, const comb_array<S>&, S**, T**, int, int, S*);
// １変数関数積プロット
template <class T, class S = double>
int plot(const pspline<T,S>&, int, S*, S**, T**, int, int = 0);
// １変数関数プロット（インデックス）
template <class T, class S = double>
int plot(const bspline<T,S>&, const comb_array<S>&, S**, T**, int, int = 0);
// １変数関数プロット
template <class T, class S = double>
int plot(const bspline<T,S>&, int, S*, S**, T**, int, int = 0);

#endif
