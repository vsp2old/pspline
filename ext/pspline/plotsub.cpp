#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include "bspline_Config.h"
#include "basis/util.h"
#include "basis/poly_array.h"
#include "basis/basis.h"
#include "basis/pspline.h"

template <typename T, typename S>
int plot(const pspline<T,S>& B, const comb_array<S>& x, S **Xp, T **Yp, int Dp, const int *b)
{
	const Int& n = x;
	int N = n.grid(), K = B.Unit_size();
	poly_index<S> H(x, Dp);
	
	int Npx = H.index_size(), L = Npx;
	for (int i = 0; i < N; i++) Xp[i] = T_ALLOC(S, Npx);
	for (int i = 0; i < K; i++) Yp[i] = T_ALLOC(T, Npx);
	while (L--) {
		int Np = H.index_value();
		poly<S> Sp = H.get_value();
		for (int i = 0; i < N; i++) Xp[i][Np] = Sp[i];
		poly_array<T> Bp = B(Sp, b);
		for (int i = 0; i < K; i++) Yp[i][Np] = Bp[i];
		++H;
	}
	return Npx;
}

template <typename T, typename S>
int plot(const pspline<T,S>& B, const comb_array<S>& x, S **Xp, T **Yp, int Dp, int Jbn, S *jb)
{
	const Int& n = x;
	int N = n.grid(), K = B.Unit_size();
	poly_index<S> H(x, Dp);
	
	int Npx = H.index_size(), L = Npx;
	for (int i = 0; i < N; i++) Xp[i] = T_ALLOC(S, Npx);
	for (int i = 0; i < K; i++) Yp[i] = T_ALLOC(T, Npx);
	while (L--) {
		int Np = H.index_value();
		poly<S> Sp = H.get_value();
		for (int i = 0; i < N; i++) Xp[i][Np] = Sp[i];
		poly_array<T> Bp = B(Sp, Jbn, jb);
		for (int i = 0; i < K; i++) Yp[i][Np] = Bp[i];
		++H;
	}
	return Npx;
}

template <typename T, typename S>
int plot(const pspline<T,S>& B, int n, S *x, S **Xp, T **Yp, int Dp, int b)
{
	int K = B.Unit_size();
	int jb[] = {b};

	int Npx = (n - 1) * Dp + 1;
	*Xp = T_ALLOC(S, Npx);
	for (int i = 0; i < K; i++) Yp[i] = T_ALLOC(T, Npx);
	for (int L = 1; L < n; L++) {
		S DT = (x[L] - x[L-1]) / Dp;
		int N0 = (L == 1 ? 0 : 1);
		for (int Nd = N0; Nd <= Dp; Nd++) {
			int Np = Nd + (L-1) * Dp;
			S Tp = x[L-1] + Nd * DT;
			poly<S> Sp(1, Tp, (S)0);
			(*Xp)[Np] = Tp;
			poly_array<T> Z = B(Sp, jb);
			for (int i = 0; i < K; ++i) Yp[i][Np] = Z[i];
		}
	}
	return Npx;
}

template <typename T, typename S>
int plot(const bspline<T,S>& B, const comb_array<S>& x, S **Xp, T **Yp, int Dp, int Jbn)
{
#ifndef NDEBUG
	const Int& n = x;
	int N = n.grid();
#endif
	assert(N == 1);
	int K = B.Unit_size();
	poly_index<S> H(x, Dp);
	
	int Npx = H.index_size(), L = Npx;
	Xp[0] = T_ALLOC(S, Npx);
	for (int i = 0; i < K; i++) Yp[i] = T_ALLOC(T, Npx);
	while (L--) {
		int Np = H.index_value();
		poly<S> Sp = H.get_value();
		(*Xp)[Np] = Sp[0];
		poly_array<T> Bp = B(Sp[0], Jbn);
		for (int i = 0; i < K; i++) Yp[i][Np] = Bp[i];
		++H;
	}
	return Npx;
}

template <typename T, typename S>
int plot(const bspline<T,S>& B, int n, S *x, S **Xp, T **Yp, int Dp, int b)
{
	int K = B.Unit_size();

	int Npx = (n - 1) * Dp + 1;
	*Xp = T_ALLOC(S, Npx);
	for (int i = 0; i < K; i++) Yp[i] = T_ALLOC(T, Npx);
	for (int L = 1; L < n; L++) {
		S DT = (x[L] - x[L-1]) / Dp;
		int N0 = (L == 1 ? 0 : 1);
		for (int Nd = N0; Nd <= Dp; Nd++) {
			int Np = Nd + (L-1) * Dp;
			S Tp = x[L-1] + Nd * DT;
			(*Xp)[Np] = Tp;
			poly_array<T> Z = B(Tp, b);
			for (int i = 0; i < K; ++i) Yp[i][Np] = Z[i];
		}
	}
	return Npx;
}

template int plot(const pspline<double,double>&, const comb_array<double>&, double**, double**, int, const int*);
template int plot(const pspline<double,double>&, const comb_array<double>&, double**, double**, int, int, double*);
template int plot(const pspline<double,double>&, int, double *, double **, double **, int, int);
template int plot(const bspline<double,double>&, const comb_array<double>&, double **, double **, int, int);
template int plot(const bspline<double,double>&, int, double *, double **, double **, int, int);

