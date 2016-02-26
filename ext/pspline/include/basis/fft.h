#ifndef _FFT_H_INCLUDED_
#define _FFT_H_INCLUDED_

template <typename T> bool rft(T *, size_t, int);	// 実数値高速フーリエ変換
template <typename T> bool cft(T *, size_t, int);	// 複素高速フーリエ変換（実数値配列引数）

template <typename T = double, typename S = T> struct fft
{

static bool complex_forward(T *dat, int N)  { return cft(dat, N << 1, 1);  }

static bool complex_backward(T *dat, int N) { return cft(dat, N << 1, 0);  }

static bool complex_inverse(T *dat, int N)  { return cft(dat, N << 1, -1); }

static bool complex_transform(T *dat, int N, int f) { return cft(dat, N << 1, f); }

static void half_en_pack(double *data, size_t n, double *dfft)
{
	size_t nn, i, j;
	double *data1 = data, *data2 = data + n;
	
	nn = n + n;
	data1[0] = dfft[0]; data2[0] = dfft[1];
	for (i = 2; i < n; i += 2) {
		j = nn - i;
		data1[i-1] = 0.5*(dfft[ j ] + dfft[ i ]);
		data1[ i ] = 0.5*(dfft[i+1] - dfft[j+1]);
		data2[i-1] = 0.5*(dfft[i+1] + dfft[j+1]);
		data2[ i ] = 0.5*(dfft[ j ] - dfft[ i ]);
	}
	if (i == n) {
		data1[i-1] = dfft[i]; data2[i-1] = dfft[i+1];
	}
}

static void complex_get(T *data, int N, int i, T *val)
{
	i %= N;
	if (N/2 < i ) i -= N;
	if (i < -N/2) i += N;
	T *var = &data[(i < 0 ? i + N : i) << 1];
	val[0] = var[0];
	val[1] = var[1];
}
	
static bspline<T,S> *cfft(T *dat, int N, int j)
{
	T *data = new T[N << 1];
	S *x = new S[N+1];
	for (int i = 0; i < N; ++i) {
		x[i] = i - N/2;
		int r = i < N/2 ? i + N/2 + N%2 : i - N/2;
		data[2*i]   = dat[r*2];
		data[2*i+1] = dat[r*2+1];
	}	x[N] = N - N/2;
	bspline<T,S> *wfft = make_bspline<2,T,S>(x, N, data, j, 1);
	delete[] x;
	delete[] data;
	return wfft;
}

static bool real_forward(T *dat, int N) { return rft(dat, N, 1); }

static bool real_backward(T *dat, int N) { return rft(dat, N, 0); }

static bool real_inverse(T *dat, int N) { return rft(dat, N, -1); }

static bool real_transform(T *dat, int N, int f) { return rft(dat, N, f); }

static void half_de_pack(double *dfft, size_t n, double *data)
{
	size_t nn, i, j;
	double *data1 = data, *data2 = data + n;
	
	nn = n + n;
	dfft[0] = data1[0]; dfft[1] = data2[0];
	for (i = 2; i < n; i += 2) {
		j = nn - i;
		dfft[ i ] = data1[i-1] - data2[ i ];
		dfft[i+1] = data2[i-1] + data1[ i ];
		dfft[ j ] = data1[i-1] + data2[ i ];
		dfft[j+1] = data2[i-1] - data1[ i ];
	}
	if (i == n) {
		dfft[n] = data1[n-1]; dfft[n+1] = data2[n-1];
	}
}

static void real_get(T *data, int N, int i, T *val)
{
	i %= N;
	if (N/2 < i ) i -= N;
	if (i < -N/2) i += N;
	int j = (i < 0 ? -i : i) << 1;
	T *var = &data[i == 0 ? 0 : j - 1];
	val[0] = var[0];
	if (i == 0 || j == N)
		val[1] = 0;
	else
		val[1] = var[1] * (i < 0 ? -1 : 1);
}

static void half_get(T *data1, int N, int i, T *val)
{
	i %= N;
	if (N/2 < i ) i -= N;
	if (i < -N/2) i += N;
	int j = (i < 0 ? -i : i) << 1;
	T *data2 = data1 + N;
	if (i == 0) {
		val[0] = data1[0];
		val[1] = data2[0];
	} else if (j == N) {
		val[0] = data1[N-1];
		val[1] = data2[N-1];
	} else if (i < 0) {
		val[0] = data1[j-1] + data2[j];
		val[1] = data2[j-1] - data1[j];
	} else {
		val[0] = data1[j-1] - data2[j];
		val[1] = data2[j-1] + data1[j];
	}
}

static bspline<T,S> *rfft(T *dat, int N, int j)
{
	T *data = new T[N << 1];
	S *x = new S[N+1];
	for (int i = 0; i < N; ++i) {
		x[i] = i - N/2;
		int r = i < N/2 ? N/2 - i : i - N/2;
		if (i == 0 && N%2 == 0) {
			data[0] = dat[2*r-1];
			data[1] = 0.0;
		} else if (i == N/2) {
			data[2*i] = dat[0];
			data[2*i+1] = 0.0;
		} else {
			data[2*i]   = dat[2*r-1];
			data[2*i+1] = dat[2*r] * (i < N/2 ? -1 : 1);
		}
	}	x[N] = N - N/2;
	bspline<T,S> *wfft = make_bspline<2,T,S>(x, N, data, j, 1);
	delete[] x;
	delete[] data;
	return wfft;
}

};	// end of fft;

#endif
