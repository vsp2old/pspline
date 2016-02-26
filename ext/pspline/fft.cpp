#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "bspline.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#define SDIM 10

template <typename T>
inline void relat(T & xr, T & xi, T wqr, T wqi)
{
	T wtemp = xr;
	xr = wtemp * wqr - xi * wqi;
	xi = wtemp * wqi + xi * wqr;
}
template <typename T>
inline void relat(T & yr, T & yi, T xr, T xi, T gr, T gi)
{
	T wtemp = yr;
	yr = wtemp * xr - yi * xi + gr;
	yi = wtemp * xi + yi * xr + gi;
}
template <typename T>
inline void relat(T *xr, T *xi, T *yr, T *yi, T wr, T wi)
{
	T tr = wr * (*yr) - wi * (*yi);
	T ti = wr * (*yi) + wi * (*yr);
	*yr  = *xr - tr; *xr += tr;
	*yi  = *xi - ti; *xi += ti;
}
/*******************************************************************************
	複素高速フーリエ変換（実数配列の引数）
	Fast Fourier Transformation / Cooley-Tukey Method
*******************************************************************************/
template <typename T>
void cft2(T *data, size_t n, int f)
{
	size_t mmax, m, j, istep, i, N = n >> 1;
	T wr, wi, wpr, wpi, theta;

	// ビット反転アルゴリズム
	for (j = 1, i = 1; i < n; i += 2) {
		if (j > i) {							// 複素数を交換
			wr = data[j-1]; data[j-1] = data[i-1]; data[i-1] = wr;
			wi = data[ j ]; data[ j ] = data[ i ]; data[ i ] = wi;
		}
		m = N;
		while(m >= 2 && j > m) { j -= m; m >>= 1; }
		j += m;
	}
	T F = -f;
	mmax = 2; theta = M_PI * F;
	while (n > mmax) {
		istep = mmax << 1;
		wpr = cos(theta); wpi = sin(theta);		// 三角関数の漸化式の初期値
		wr = 1.0; wi = 0.0;
		for (m = 1; m < mmax; m += 2) { 		// 2重の内側ループ
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				relat<T>(&data[i-1], &data[i], &data[j-1], &data[j], wr, wi);
			}
			relat<T>(wr, wi, wpr, wpi); 		// 三角関数の漸化式
		}
		mmax = istep;
		theta /= 2.0;
	}
}
/*
	任意基数のFFT
*/
template <typename T, size_t R>
void cft(T *data, size_t n, int f)
{
	size_t h, i, j, k, L, m, mmax, istep;
	size_t N = n >> 1, Nr = N / R;
	T wr, wi, wpr, wpi, wqr, wqi, Xr, Xi, Yr, Yi, *gw;
	gw = new T[R*2];
	// R進数の桁反転
	for (j = 1, i = 1; i < n; i += 2) {
		if (i < j) {
			wr = data[i-1]; data[i-1] = data[j-1]; data[j-1] = wr;
			wi = data[ i ]; data[ i ] = data[ j ]; data[ j ] = wi;
		}
		m = (R-1) * (Nr << 1);
		while (m >= 2 && j > m) {j -= m; m /= R;}
		j += m / (R-1);
	}
	T F = -f;
	T theta = 2 * M_PI * F / (T)R;
	wqr = cos(theta); wqi = sin(theta);
	mmax = 2;
	while (n > mmax) {
		istep = mmax * R;
		wpr = cos(theta); wpi = sin(theta);
		wr = 1.0; wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				for (j = i, L = 0; L < R; L++, j += mmax) {
					h = (L << 1) + 1;
					gw[h-1] = data[j-1]; gw[h] = data[j];
				}
				// バタフライ演算
				Xr = wr; Xi = wi;
				for (j = i, L = 0; L < R; L++, j += mmax) {
					Yr = Yi = 0;
					for (k = R; k > 0; --k) {
						h = (k << 1) - 1;
						relat<T>(Yr, Yi, Xr, Xi, gw[h-1], gw[h]);
					}	data[j-1] = Yr; data[j] = Yi;
					relat<T>(Xr, Xi, wqr, wqi);
				}
			}
			relat<T>(wr, wi, wpr, wpi);
		}
		mmax = istep;
		theta /= R;
	}
	delete[] gw;
}
/*
	複数基数のFFT subroutine
*/
template <typename T>
void cft(size_t r, T *data, size_t n, int f, size_t stride = 1)
{
	size_t h, i, j, k, L, Li, Lj;
	T xr, xi, yr, yi, *gw;
	size_t N = n >> 1, Nr = N / r;	// N = r**M;
	gw = new T[r*2];
	// r進数の桁反転
	for (j = Nr, i = 1; i < N - 1; ++i) {
		if (i < j) {
			size_t ki = (i * stride) << 1;
			size_t kj = (j * stride) << 1;
			xr = data[ ki ]; data[ ki ] = data[ kj ]; data[ kj ] = xr;
			xi = data[ki+1]; data[ki+1] = data[kj+1]; data[kj+1] = xi;
		}
		size_t m = (r-1) * Nr;
		while (m > 0 && j >= m) {j -= m; m /= r;}
		j += m / (r-1);
	}
	T F = -f;
	T theta = 2 * M_PI * F / r;
	T wqr = cos(theta), wqi = sin(theta);
	size_t Nj = 1;
	while (Nj < N) {
		size_t step = Nj * r;
		T wpr = cos(theta), wpi = sin(theta);
		T wr = 1.0, wi = 0.0;
		for (L = 0; L < Nj; ++L) {
			for (h = L; h < N; h += step) {
				size_t m = h * stride;
				for (Li = m, Lj = 0; Lj < r; ++Lj, Li += Nj * stride) {
					i = Li << 1; j = Lj << 1;
					gw[j] = data[i]; gw[j+1] = data[i+1];
				}
				// バタフライ演算
				xr = wr; xi = wi;
				for (Li = m, Lj = 0; Lj < r; Lj++, Li += Nj * stride) {
					yr = yi = 0;
					for (k = r; k > 0; --k) {
						j = (k << 1) - 1;
						relat<T>(yr, yi, xr, xi, gw[j-1], gw[j]);
					}	i = Li << 1;
						data[i] = yr; data[i+1] = yi;
					relat<T>(xr, xi, wqr, wqi);
				}
			}
			relat<T>(wr, wi, wpr, wpi);
		}
		theta /= r;
		Nj = step;
	}
	delete[] gw;
}
/*
	複数基数のFFT main routine
	g : データ配列
	N : データ個数
	f : 1 forward, 0 backword, -1 invert
*/
template <typename T>
bool cft(T *data, size_t n, int f)
{
	size_t N = n >> 1, Nres = N, s = 0, k;
	size_t Ni[SDIM], ni[SDIM], M[SDIM], r[SDIM] = {2};
	// Chinese Remainder Theorem
	do {
		M[s] = 0; ni[s] = 1; Ni[s] = N;
		while (Ni[s] % r[s] == 0) {
			 M[s]++;
			ni[s] *= r[s];	// ni = r**M;
			Ni[s] /= r[s];
		}	Nres /= ni[s];
		int rnew;
		if (Nres != 1) {
			rnew = r[s] + (r[s] == 2 ? 1 : 2);
			if (M[s] > 0) s++;
			if (s < SDIM) {
				// 素因数分解
				while (Nres % rnew != 0) rnew += 2;
				r[s] = rnew;
			}
		}
	} while (Nres > 1 && s < SDIM);
	int F = (f > 0) ? 1 : -1;
	if (s >= SDIM) return false; else
	if (s == 0) {
		switch (r[s]) {
			case 2: cft2<T>(data, n, F); break;
			case 3: cft<T,3>(data, n, F); break;
			case 5: cft<T,5>(data, n, F); break;
			case 7: cft<T,7>(data, n, F); break;
			default: cft<T>(r[s], data, n, F);
		}
	} else {
		size_t smax = s + 1;
		T gw[n];
		// 最初の並び替え
		for (k = 0; k < N; ++k) {
			size_t ks, kk = 0;
			for (s = 0; s < smax; ++s)
				kk = kk * ni[s] + k % ni[s];
				kk <<= 1; ks = k << 1;
			 gw[kk] = data[ks]; gw[kk+1] = data[ks+1];
		}
		size_t step, stride = 1;
		// フーリエ変換
		for (k = smax; k > 0; --k) {
			s = k - 1;
			step = stride * ni[s];
			for (size_t k1 = 0; k1 < stride; ++k1)
				for (size_t k0 = k1; k0 < N; k0 += step)
					cft(r[s], &gw[k0<<1], ni[s] << 1, F, stride);
			stride = step;
		}
		// 最終の並び替え
		for (k = 0; k < N; ++k) {
			size_t Ls, Lu, L = 0, LL = 0;
			for (s = 0; s < smax; ++s) {
				Ls = k % ni[s];
				LL = LL * ni[s] + Ls;
				L += Ni[s] * Ls;
			}	Ls = (L % N) << 1;
				Lu = LL << 1;
			data[Ls] = gw[Lu]; data[Ls+1] = gw[Lu+1];
		}
	}
	if (f < 0) for (k = 0; k < n; ++k) data[k] /= N;
	return true;
}
/*******************************************************************************
	実数値高速フーリエ変換
	Fast Fourier Transformation / Cooley-Tukey Method
*******************************************************************************/
/*
	基数２のFFT
*/
template <typename T>
void rft2(T *data, size_t N, int f)
{
	size_t i, j, NH = N >> 1;
	T c, g1r, g1i, g2r, g2i, temp;
	T wr, wi, wpr, wpi, wtemp;
	T F = -f, theta = M_PI * F / (T)NH;

	wpr = cos(theta); wpi = sin(theta);
	if (f < 0) {
		wr = data[0]; wi = data[N-1];
		data[0] = wr + wi; wtemp = wr - wi;
		wr = 0.0; wi = 1.0;
		for (i = 2; i <= NH; i += 2) {
			relat<T>(wr, wi, wpr, wpi);
			j = N - i; temp = data[i-1]; data[i-1] = wtemp;
			g1r = (temp + data[j-1]); g1i = (data[i] - data[j]);
			g2r = (temp - data[j-1]); g2i = (data[i] + data[j]);
			relat<T>(g2r, g2i, wr, wi);
			if (i < j) {
				data[i] = g1r + g2r; wtemp = g1i + g2i;
			}
			data[j] = g1r - g2r; data[j+1] = g2i - g1i;
		}
	}
	cft2(data, N, f);
	if (f > 0) {
		c = 0.5;
		wr = data[0]; wi = data[1];
		data[0] = wr + wi; wtemp = wr - wi;
		wr = 0.0; wi = -1.0;
		for (i = 2; i <= NH; i += 2) {
			relat<T>(wr, wi, wpr, wpi);
			j = N - i; temp = data[j+1]; data[j+1] = wtemp;
			g1r = c*(data[i] + data[j]); g1i = c*(data[i+1] - temp);
			g2r = c*(data[i] - data[j]); g2i = c*(data[i+1] + temp);
			relat<T>(g2r, g2i, wr, wi);
			data[i-1] = g1r + g2r; data[i] = g1i + g2i;
			if (i < j) {
				wtemp = g1r - g2r; data[j] = g2i - g1i;
			}	
		}
	}
}
/*
	任意基数のFFT
*/
template <typename T, size_t R>
void rft(T *data, size_t N, int f)
{
	size_t h, i, j, k, L, m, mmax, istep;
	size_t n = N << 1, Nr = N / R;
	T wr, wi, wpr, wpi, wqr, wqi, *temp;
	T xr, xi, yr, yi, *gw;
	temp = new T[n]; gw = new T[R*2];
	// R進数の桁反転
	for (j = 1, i = 0; i < n; i += 2) {
		if (f > 0) {
			xr = data[i>>1];
			xi = 0;
		} else if (i == 0) {
			xr = data[0];
			xi = 0;
		} else if (i < n/2) {
			xr = data[i-1];
			xi = data[ i ];
		} else if (i == n/2) {
			xr = data[i-1];
			xi = 0;
		} else {
			xr = data[n-i-1];
			xi =-data[ n-i ];
		}
		temp[j-1] = xr; temp[j] = xi;
		m = (R-1) * (Nr << 1);
		while (m >= 2 && j > m) {j -= m; m /= R;}
		j += m / (R-1);
	}
	T F = -f;
	T theta = 2.0 * M_PI * F / (T)R;
	wqr = cos(theta); wqi = sin(theta);
	mmax = 2;
	while (n > mmax) {
		istep = mmax * R;
		wpr = cos(theta); wpi = sin(theta);
		wr = 1.0; wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				for (j = i, L = 0; L < R; L++, j += mmax) {
					h = (L << 1) + 1;
					gw[h-1] = temp[j-1]; gw[h] = temp[j];
				}
				// バタフライ演算
				xr = wr; xi = wi;
				for (j = i, L = 0; L < R; L++, j += mmax) {
					yr = yi = 0;
					for (k = R; k > 0; --k) {
						h = (k << 1) - 1;
						relat<T>(yr, yi, xr, xi, gw[h-1], gw[h]);
					}	temp[j-1] = yr; temp[j] = yi;
					relat<T>(xr, xi, wqr, wqi);
				}
			}
			relat<T>(wr, wi, wpr, wpi);
		}
		mmax = istep;
		theta /= R;
	}
	data[0] = temp[0];
	for (i = 1; i < N; ++i) data[i] = temp[f > 0 ? i + 1 : i << 1];
	delete[] gw; delete[] temp;
}
/*
	複数基数のFFT subroutine
*/
template <typename T>
void rft(size_t r, T *data, size_t N, int f)
{
	size_t h, i, j, k, L, m, mmax, istep;
	size_t n = N << 1, Nr = N / r;
	T wr, wi, wpr, wpi, wqr, wqi, *temp;
	T xr, xi, yr, yi, *gw;
	temp = new T[n]; gw = new T[r*2];
	// r進数の桁反転
	for (j = 1, i = 0; i < n; i += 2) {
		if (f > 0) {
			xr = data[i>>1];
			xi = 0;
		} else if (i == 0) {
			xr = data[0];
			xi = 0;
		} else if (i < N) {
			xr = data[i-1];
			xi = data[ i ];
		} else if (i == N) {
			xr = data[i-1];
			xi = 0;
		} else {
			xr = data[n-i-1];
			xi =-data[ n-i ];
		}
		temp[j-1] = xr; temp[j] = xi;
		m = (r-1) * (Nr << 1);
		while (m >= 2 && j > m) {j -= m; m /= r;}
		j += m / (r-1);
	}
	T F = -f;
	T theta = 2.0 * M_PI * F / (T)r;
	wqr = cos(theta); wqi = sin(theta);
	mmax = 2;
	while (n > mmax) {
		istep = mmax * r;
		wpr = cos(theta); wpi = sin(theta);
		wr = 1.0; wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				for (j = i, L = 0; L < r; L++, j += mmax) {
					h = (L << 1) + 1;
					gw[h-1] = temp[j-1]; gw[h] = temp[j];
				}
				// バタフライ演算
				xr = wr; xi = wi;
				for (j = i, L = 0; L < r; L++, j += mmax) {
					yr = yi = 0;
					for (k = r; k > 0; --k) {
						h = (k << 1) - 1;
						relat<T>(yr, yi, xr, xi, gw[h-1], gw[h]);
					}	temp[j-1] = yr; temp[j] = yi;
					relat<T>(xr, xi, wqr, wqi);
				}
			}
			relat<T>(wr, wi, wpr, wpi);
		}
		mmax = istep;
		theta /= r;
	}
	data[0] = temp[0];
	for (i = 1; i < N; ++i) data[i] = temp[f > 0 ? i + 1 : i << 1];
	delete[] gw; delete[] temp;
}
/*
	複数基数のFFT main routine
	g : データ配列
	N : データ個数
	f : 1 forward, 0 backword, -1 invert
*/
template <typename T>
bool rft(T *data, size_t N, int f)
{
	size_t Nres = N, s = 0, k, n = N << 1;
	size_t Ni[SDIM], ni[SDIM], M[SDIM], r[SDIM]; r[0] = 2;
	// Chinese Remainder Theorem
	do {
		M[s] = 0; ni[s] = 1; Ni[s] = N;
		while (Ni[s] % r[s] == 0) {
			 M[s]++;
			ni[s] *= r[s];	// ni = r**M;
			Ni[s] /= r[s];
		}	Nres /= ni[s];
		int rnew;
		if (Nres != 1) {
			rnew = r[s] + (r[s] == 2 ? 1 : 2);
			if (M[s] > 0) s++;
			if (s < SDIM) {
				// 素因数分解
				while (Nres % rnew != 0) rnew += 2;
				r[s] = rnew;
			}
		}
	} while (Nres > 1 && s < SDIM);
	int F = (f > 0) ? 1 : -1;
	if (s >= SDIM) return false; else
	if (s == 0) {
		switch (r[s]) {
			case 2: rft2<T>(data, N, F); break;
			case 3: rft<T,3>(data, N, F); break;
			case 5: rft<T,5>(data, N, F); break;
			case 7: rft<T,7>(data, N, F); break;
			default: rft<T>(r[s], data, N, F);
		}
	} else {
		size_t smax = s + 1;
		T gw[n];
		// 最初の並び替え
		for (k = 0; k < N; ++k) {
			size_t ks, kk = 0;
			T wr, wi;
			for (s = 0; s < smax; ++s)
				kk = kk * ni[s] + k % ni[s];
				kk <<= 1; ks = k << 1;
			if (F > 0) {
				wr = data[k];
				wi = 0;
			} else if (ks == 0) {
				wr = data[0];
				wi = 0;
			} else if (ks == N) {
				wr = data[N-1];
				wi = 0;
			} else if (ks < N) {
				wr = data[ks-1];
				wi = data[ ks ];
			} else {
				ks = (N-k) << 1;
				wr = data[ks-1];
				wi =-data[ ks ];
			}
			gw[kk] = wr; gw[kk+1] = wi;
		}
		size_t step, stride = 1;
		// フーリエ変換
		for (k = smax; k > 0; --k) {
			s = k - 1;
			step = stride * ni[s];
			for (size_t k1 = 0; k1 < stride; ++k1)
				for (size_t k0 = k1; k0 < N; k0 += step)
					cft(r[s], &gw[k0<<1], ni[s] << 1, F, stride);
			stride = step;
		}
		// 最終の並び替え
		for (k = 0; k < N; ++k) {
			size_t Ls, Lu, L = 0, LL = 0;
			for (s = 0; s < smax; ++s) {
				Ls = k % ni[s];
				LL = LL * ni[s] + Ls;
				L += Ni[s] * Ls;
			}	Ls = (L % N);
				Lu = LL << 1;
			if (F > 0) {
				Ls <<= 1;
				if (Ls == 0) data[Ls] = gw[Lu];
				else
				if (Ls == N) data[Ls-1] = gw[Lu];
				else
				if (Ls < N) {
					data[Ls-1] = gw[ Lu ];
					data[ Ls ] = gw[Lu+1];
				}
			} else  data[ Ls ] = gw[ Lu ];
		}
	}
	if (f < 0) for (k = 0; k < N; ++k) data[k] /= N;
	return true;
}

template bool cft<double>(double *, size_t, int);
template bool rft<double>(double *, size_t, int);
