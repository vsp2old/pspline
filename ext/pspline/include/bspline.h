#ifndef _BSPLINE_H_INCLUDED_
#define _BSPLINE_H_INCLUDED_

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "bspline_Config.h"
#include "basis/util.h"
#include "basis/poly_array.h"
#include "basis/basis.h"
#include "basis/pspline.h"
#include "basis/uspline.h"

/*******************************************************************************
	make object utility

	1. poly<T,N> = <{TN,...,T1},T0>

*******************************************************************************/

// make_alias<T,N=1>(z) = @<{z0,...,zN-1},zN>

template<typename T, int N = 1> inline
poly<T> make_alias(T *z)
{
	return poly<T>(z, N);
}

// make_atom<T>(a) = @<a0>

template <typename T>
poly<T> make_atom(T *a)
{
	return poly<T>(a, 0);
}

// make_poly<T,tN,...,t1>(t0=1) = <{tN,...,t1},t0>

template<typename T, T... args> inline
poly<T> make_poly(T a = 1)
{
	int N = sizeof...(args);
	T s[] = {args...};
	if (N == 0)
	return poly<T>(0, a);
	else
	return poly<T>(N, s, a);
}

// make_poly<T,N>(tN,...,t1,t0) = <{tN,...,t1},t0>

template<typename T, int N, class... Args> inline
poly<T> make_poly(T car, Args... cdr)
{
	return poly<T>(N, car, cdr...);
}

// make_poly<T,N>(s) = <{s0,...,sN-1},sN>

template<typename T, int N = 1> inline
poly<T> make_poly(const T *s)
{
	return poly<T>(N, s);
}

// make_poly<T,N>(s, a0) = <{s0,...,sN-1},a0>

template<typename T, int N = 1> inline
poly<T> make_poly(const T *s, T a)
{
	return poly<T>(N, s, a);
}

/*******************************************************************************

	2. poly_array<T,N> = {<{SN,...,S1},S0>,[T...]}

*******************************************************************************/

// make_array<T,SN,...,S1,S0>(t) = [<{SN,...,S1},S0>,[t0,...]]

template<class T, int... args> inline
poly_array<T> make_array(const T *dat)
{
	int N = sizeof...(args) - 1;
	int s[] = {args...};
	return poly_array<T>(dat, N, s);
}

// make_array<T,N>(t,SN,...,S1,S0) = [<{SN,...,S1},S0>,[t0,...]]

template<class T, int N, class... Args> inline
poly_array<T> make_array(const T *dat, int car, Args... cdr)
{
	int i = sizeof...(cdr), n = i - N;
	int s[] = {car, cdr...}, z[N+1];
	for (i = 0; i <= N; ++i, ++n)
		if (n < 0) z[i] = 1;
		else z[i] = s[n];
	return poly_array<T>(dat, N, z);
}

// make_array<T,N>(t,<{SN,...,S1},S0>) = [<{SN,...,S1},S0>,[t0,...]]

template<class T, int N> inline
poly_array<T> make_array(const T *dat, const poly<int>& e)
{
	return poly_array<T>(dat, e);
}

// sub_array(int i,[<{SN+1,SN,...,S1},S0>,[[t0,...],...]]) = [<{SN,...,S1},S0>,[ti,...]]

template<typename T> inline
poly_array<T> sub_array(int p, const poly_array<T>& a)
{
	return poly_array<T>(p, a);
}

// unit_array<T,N>(a) = [<N>,[a0,...,aN-1]] : N > 0
// unit_array<T>(a)   = [<1>,[a0]]          : N = 0

template<class T, int N = 0> inline
poly_array<T> unit_array(T *a)
{
	if (N > 0) return poly_array<T>(a, N);
	else	   return poly_array<T>(a);
}

// make_unit<T>(a,N) = [<N>,[a0,...,aN-1]] : N > 0
// make_unit<T>(a)   = [<1>,[a0]]          : N = 0

template<class T> inline
poly_array<T> make_unit(T *w, int N = 0)
{
	if (N > 0) return poly_array<T>(w, N);
	else	   return poly_array<T>(w);
}

/*******************************************************************************

	3. poly_view<T,N> = {<{SN,...,S1},S0>,@[T...]}

*******************************************************************************/

// make_view<T,SN,...,S1,S0>(t) = [<{SN,...,S1},S0>,@[t0,...]]

template<class T, int... args> inline
poly_view<T> make_view(T *dat)
{
	int N = sizeof...(args) - 1;
	int s[] = {args...};
	return poly_view<T>(dat, N, s);
}

// make_view<T,N>(t,SN,...,S1,S0) = [<{SN,...,S1},S0>,@[t0,...]]

template<class T, int N, class... Args> inline
poly_view<T> make_view(T *dat, int car, Args... cdr)
{
	int i = sizeof...(cdr), n = i - N;
	int s[] = {car, cdr...}, z[N+1];
	for (i = 0; i <= N; ++i, ++n)
		if (n < 0) z[i] = 1;
		else z[i] = s[n];
	return poly_view<T>(dat, N, z);
}

// make_view<T,N>(t,<{SN,...,S1},S0>) = [<{SN,...,S1},S0>,@[t0,...]]

template<class T, int N> inline
poly_view<T> make_view(T *dat, const poly<int>& e)
{
	return poly_view<T>(dat, e);
}

// view([<{SN,...,S1},S0>,[t0,...]]) = [<{SN,...,S1},S0>,@[t0,...]]

template<typename T> inline
poly_view<T> view(const poly_array<T>& c)
{
	return poly_view<T>(c);
}

// sub_view(int i,[<{SN,...,S1},S0>,[[t0,...],...]]) = [<{SN-1,...,S1},S0>,@[ti,...]]

template<typename T> inline
poly_view<T> sub_view(int p, const poly_array<T>& a)
{
	return poly_view<T>(p, a);
}

// unit_view<T,N>(a) = [<N>,@[a0,...,aN-1]] : N > 0
// unit_view<T>(a)   = [<1>,@[a0]]          : N = 0

template<class T, int N = 0> inline
poly_view<T> unit_view(T *a)
{
	if (N > 0)
	return poly_view<T>(a, N);
	else
	return poly_view<T>(a);
}

// view<T>(a,N) = [<N>,@[a0,...,aN-1]] : N > 0
// view<T>(a)   = [<1>,@[a0]]          : N = 0

template<class T> inline
poly_view<T> view(T *w, int N = 0)
{
	if (N > 0)
	return poly_view<T>(w, N);
	else
	return poly_view<T>(w);
}

/*******************************************************************************

	4. comb_array<T,N> = {<{SN,...,S1},S0>,[TN...],...,[T1,...]}

*******************************************************************************/

template<class T, int... args> inline
comb_array<T> make_comb(const T *dat)
{
	int n = sizeof...(args);
	int s[] = {args...};
	return comb_array<T>(dat, n-1, s);
}

template<class T, int N, class... args> inline
comb_array<T> make_comb(const T *dat, int b, args... arg)
{
	int s[] = {b, arg...};
	return comb_array<T>(dat, N, s);
}

template<class T, int N> inline
comb_array<T> make_comb(const T *dat, int *s)
{
	return comb_array<T>(dat, N, s);
}

template<class T, int N> inline
comb_array<T> make_comb(const T *dat, int *s, int k)
{
	return comb_array<T>(dat, N, s, k);
}

template<class T, int N, class E> inline
comb_array<T> make_comb(const T *dat, const E& e)
{
	return comb_array<T>(dat, e);
}

template<class T> inline
comb_array<T> atom_comb(const T *dat)
{
	return comb_array<T>(dat);
}

template<class T> inline
comb_array<T> atom_comb(const T *dat, int s)
{
	return comb_array<T>(dat, Atom, s);
}

/*******************************************************************************

	5. comb_view<T,N> = {<{SN,...,S1},S0>,@[T...]}

*******************************************************************************/

template<class T> inline
comb_view<T> comb(const comb_array<T>& ca)
{
	return comb_view<T>(ca);
}

template<class T, int N, class... args> inline
comb_view<T> comb(T *dat, int b, args... arg)
{
	int s[] = {b, arg...};
	return comb_view<T>(dat, N, s);
}

template<class T, int N> inline
comb_view<T> comb(T *dat, int *s)
{
	return comb_view<T>(dat, N, s);
}

template<class T, int N> inline
comb_view<T> comb(T *dat, int *s, int a)
{
	return comb_view<T>(dat, N, s, a);
}

/*******************************************************************************

	6. comp_array<T,N> = {<{ON,...,O1},{SN,...,S1},K>,[T...]}

*******************************************************************************/

template<class T>
comp_array<T> make_mono(T *dat, int o, int s, int k)
{
	return comp_array<T>(dat, o, s, k);
}

template<class T, int N, class... args>
comp_array<T> make_comp(T *dat, int b, args... arg)
{
	int s[] = {b, arg...};
	return comp_array<T>(dat, N, s);
}

template<class T, int N>
comp_array<T> make_comp(T *dat, int *s)
{
	return comp_array<T>(dat, N, s);
}

template <typename T> inline
poly_index<T> make_index(const comb_array<T>& x, int dp)
{
	return poly_index<T>(x, dp);
}

template <int M = 1, typename S = double, typename T = double> inline
bspline<T,S> *make_bspline(S *x, int n, T *y, int j, int s = 0, int c = 0, int *d = NULL)
{
	auto Y = make_array<T,1>(y, n+c, M);
	return new bspline<T,S>(Y, n+(s%2), x, j, s, d);
}

template <int N, typename T = double> inline
comb_array<T> pack_comb(const T *u, const int *n, const int *s, int a = 1)
{
	auto nc = poly<int>(N, n, a);
	auto sc = poly<int>(N, s, 0);
	auto ns = nc + sc%2;
	int nn = ns.comb_size();
	T U[nn], *p = U;
	for (int i = 0; i < N; ++i) {
		int size = ns[i];
		for (int j = 0; j < size; ++j)
			*p++ = u ? *u++ : (T)j;
	}
	return comb_array<T>(U, ns);
}

template <typename T = double> inline
comb_array<T> pack_comb(const T *u, int n, int s = 0, int a = 1)
{
	auto ns = poly<int>(1, n+s%2, a);
	int nn = ns.comb_size();
	T U[nn];
	for (int i = 0; i < nn; ++i) U[i] = u ? u[i] : (T)i;
	return comb_array<T>(U, ns);
}

template <int N, int M = 1, typename S = double, typename T = double> inline
pspline<T,S> *make_pspline(const S *x, const int *p, const T *y, const int *n, const int *j, const int *s)
{
	auto pc = poly<int>(N, p, 1);
	auto nc = poly<int>(N, n, M);
	auto jc = poly<int>(N, j, 0);
	auto sc = poly<int>(N, s, 0);
	auto U = make_comb<S,N>(x, pc);
	auto W = make_array<T,N>(y, nc);
	return new pspline<T,S>(W, U, jc, sc);
}

#include "basis/fft.h"

#endif
