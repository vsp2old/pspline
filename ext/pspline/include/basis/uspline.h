#ifndef _USPLINE_5_H_INCLUDED_
#define _USPLINE_5_H_INCLUDED_

template <typename T = double> struct Func { int jbn, jis; T val; };
template <typename T = double> struct Bound { int pos, jbn; T val; };
template <typename T = double> using Factor = T(*)(T);

/*******************************************************************************
	uspline<T>{base_spline*[N]}
*******************************************************************************/
template <typename T = double> class uspline : public base_spline<T>
{

T *data;

//	汎関数線型項
void hankansu1(const Func<T>& F, T *hr)
{
	T  *q = base_spline<T>::knots;
	int K = base_spline<T>::rank, J = K - 1;
	int Imax = base_spline<T>::imax;
	T t = q[Imax];
	T B[K];
	int kset = base_spline<T>::basic(t, K, B);
	int ks = kset - J;
	for (int l = 0; l < Imax; ++l) {
		T aj[Imax], ac[Imax];
		for (int i = 0; i < Imax; ++i) ac[i] = 0;
		ac[l] = 1;
		aj[0] = ac[0] * (q[K] - q[0]) / K;
		for (int j = 1; j < Imax; ++j)
			aj[j] = aj[j-1] + ac[j] * (q[j+K] - q[j]) / K;
		T rn = 0;
		for (int i = ks; i <= kset; ++i)
			rn += aj[i] * B[i - ks];
		hr[l] += F.val * rn;
	}
}

//	汎関数２乗項
void hankansu2(const Func<T>& F, T *hl)
{
	T  *q = base_spline<T>::knots;
	int K = base_spline<T>::rank, J = K - 1;
	int Imax = base_spline<T>::imax;
	T ld[K*2];
	for (int l = 0; l < K*2; ++l) {
		if (l == 0)
			ld[l] = 1.0;
		else
			ld[l] = ld[l-1] * l;
	}
	int jb = F.jbn;
	for (int l = 0; l < Imax; ++l) {
		int ql = l + J;
		T t = q[ql];
		parray<T> fl = base_spline<T>::bases(t, J);
		int ks = fl(0), kt = fl(1), kset = ks + J;
		T *ft = *fl;
		for (int i = ks; i <= kset; ++i) {
			for (int j = ks; j <= kset; ++j) {
				for (int r = jb; r < K; ++r) {
					int rb = r - jb;
					T Bik = ft[i - ks + r*kt];
					for (int s = jb; s < K; ++s) {
						int sb = s - jb;
						T Bjk = ft[j - ks + s*kt];
						int lm = rb + sb + 1;
						T Ql = pow(q[ql+1] - q[ql], lm);
						T p = Bik * Bjk * Ql / (ld[rb] * ld[sb] * lm);
						hl[i*Imax + j] += F.val * p;
					}
				}
			}
		}
	}
}

void gyoretu_ritz(double *hl, double *hr, double *Bd, double *Dy)
{
	int Imax = base_spline<T>::imax;
	for (int l = 0; l < Imax; ++l) {
		for (int i = 0; i < Imax; ++i)
			Bd[l*Imax + i] = 2 * hl[l*Imax + i];
			Dy[l] = -hr[l];
	}
}

void zyoken(int lcc, int *ibd, const Bound<T> *Bs, T *Bd, T *Dx, T *Dy)
{
	int K = base_spline<T>::rank, J = K - 1;
	int Imax = base_spline<T>::imax;
	T B[K];
	for (int l = 0; l < lcc; l++) {
		T t = Dx[Bs[l].pos];
		int jb = Bs[l].jbn;
		int kset = base_spline<T>::basic(t, K, B, -jb);
		int ks = kset - J;
		for (int i = 0; i < Imax; i++)
			Bd[ibd[l]*Imax + i] = (i < ks || kset < i) ? 0 : B[i - ks];
			Dy[ibd[l]] = Bs[l].val;
	}
}

public:
uspline() : data(NULL) {}
~uspline() {delete[] data;}
uspline(int R, Func<T> *df, int n, T *x, int j, int d, int lcc, Bound<T> *bf, int *lbd)
 : base_spline<T>(n, x, j, d)
{
	int Imax = base_spline<T>::imax;
	T hl[Imax * Imax], hr[Imax];
	T Bd[Imax][Imax], Dy[Imax];
	for (int i = 0; i < Imax; ++i) {
		for (int j = 0; j < Imax; ++j)
			hl[i*Imax+j] = 0;
			hr[i] = 0;
	}
	for (int term = 0; term < R; ++term)
		if (df[term].jis == 2)
			hankansu2(df[term], hl);
		else
			hankansu1(df[term], hr);
	gyoretu_ritz(hl, hr, (T*)Bd, Dy);
	zyoken(lcc, lbd, bf, (T*)Bd, x, Dy);
	data = new T[Imax];
	for (int i = 0; i < Imax; ++i) data[i] = Dy[i];
	lu_solve((T*)Bd, Imax, data);
}

T operator [] (T t) const { return base_spline<T>::apply(t, data); }

}; // end of uspline;

/*******************************************************************************
	nspline<T>{base_spline*[N]}
*******************************************************************************/
template <typename T> class nspline : public base_spline<T>
{
T *data;

//	微分方程式左辺
void sa_hen(T *Dx, int term, Func<T> *F, Factor<T> *C, T *Bd)
{
	int K = base_spline<T>::rank, J = K - 1;
	int Imax = base_spline<T>::imax;
	T B[K*(K+1)];
	for (int l = 0; l < Imax; ++l) {
		T t = Dx[l];
		for (int k = 0; k < term; ++k) {
			F[k].val = C[k](t);
			int jb = F[k].jbn;
			int kset = base_spline<T>::basic(t, K, B, -jb);
			int ks = kset - J;
			for (int i = 0; i < Imax; ++i)
				Bd[l*Imax + i] += F[k].val * ((ks <= i && i <= kset) ? B[i - ks] : 0.0);
		}
	}
}

//	微分方程式右辺
void u_hen(Factor<T>& f, int Imax, T *Dx, T *Dy)
{
	for (int i = 0; i < Imax; ++i) Dy[i] = f(Dx[i]);
}

void zyoken(int lcc, int *ibd, const Bound<T> *Bs, T *Bd, T *Dx, T *Dy)
{
	int K = base_spline<T>::rank, J = K - 1;
	int Imax = base_spline<T>::imax;
	T B[K];
	for (int l = 0; l < lcc; l++) {
		T t = Dx[Bs[l].pos];
		int jb = Bs[l].jbn;
		int kset = base_spline<T>::basic(t, K, B, -jb);
		int ks = kset - J;
		for (int i = 0; i < Imax; i++)
			Bd[ibd[l]*Imax + i] = (i < ks || kset < i) ? 0 : B[i - ks];
			Dy[ibd[l]] = Bs[l].val;
	}
}

public:
nspline() : data(NULL) {}
~nspline(){delete[] data;}
nspline(int R, Func<T> *df, Factor<T> *uf, int n, T *x, int j, int d, int lcc, Bound<T> *bf, int *lbd)
 : base_spline<T>(n, x, j, d)
{
	int Imax = base_spline<T>::imax;
	double Bd[Imax * Imax], Dy[Imax];
	for (int i = 0; i < Imax; ++i) {
		for (int j = 0; j < Imax; ++j)
			Bd[i*Imax+j] = 0;
			Dy[i] = 0;
	}
	sa_hen(x, R, df, uf, Bd);
	u_hen(uf[R], Imax, x, Dy);
	zyoken(lcc, lbd, bf, Bd, x, Dy);
	data = new T[Imax];
	for (int i = 0; i < Imax; ++i) data[i] = Dy[i];
	lu_solve(Bd, Imax, data);
}

T operator [] (T t) const { return base_spline<T>::apply(t, data); }

T operator () (T t, int Jbn = 0) const { return base_spline<T>::apply(t, data, Jbn); }

}; // end of nspline

#endif
