/*
Copyright 2025, Yves Gallot

mersenne2.cpp is free source code. You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>
#include <ctime>

// Z/{2^61 - 1}Z: the prime field of order p = 2^61 - 1
class Z61
{
private:
	uint64_t _n;

	static const uint64_t _p = (uint64_t(1) << 61) - 1;
	static const uint64_t _primroot = 37;

public:
	Z61() {}
	explicit Z61(const uint64_t n) : _n(n) {}

	uint64_t get() const { return _n; }

	bool operator==(const Z61 & rhs) const { return (_n == rhs._n); }

	Z61 operator-() const { return Z61((_n == 0) ? 0 : _p - _n); }

	Z61 & operator+=(const Z61 & rhs) { const uint64_t c = (_n >= _p - rhs._n) ? _p : 0; _n += rhs._n; _n -= c; return *this; }
	Z61 & operator-=(const Z61 & rhs) { const uint64_t c = (_n < rhs._n) ? _p : 0; _n -= rhs._n; _n += c; return *this; }
	Z61 & operator*=(const Z61 & rhs) { _n = uint64_t((_n * __uint128_t(rhs._n)) % _p); return *this; }

	Z61 operator+(const Z61 & rhs) const { Z61 r = *this; r += rhs; return r; }
	Z61 operator-(const Z61 & rhs) const { Z61 r = *this; r -= rhs; return r; }
	Z61 operator*(const Z61 & rhs) const { Z61 r = *this; r *= rhs; return r; }

	Z61 pow(const size_t e) const
	{
		if (e == 0) return Z61(1u);

		Z61 r = Z61(1u), y = *this;
		for (size_t i = e; i != 1; i /= 2)
		{
			if (i % 2 != 0) r *= y;
			y *= y;
		}
		r *= y;

		return r;
	}

	Z61 invert() const { return pow(_p - 2); }

	static const Z61 root_nth(const size_t n) { return Z61(_primroot).pow((_p - 1) / n); }
};

static const Z61 J = Z61::root_nth(3), inv2 = Z61(2).invert(), inv3 = Z61(3).invert();

// GF((2^61 - 1)^2): the prime field of order p^2, p = 2^61 - 1
class GF61
{
private:
	Z61 _s0, _s1;
	// a primitive root of order 3*2^62 which is a root of (0, 1).
	static const uint64_t _h_order = uint64_t(3) << 62;
	static const uint64_t _h_0 = 1794883824738941ull, _h_1 = 7671273768663199ull;

public:
	GF61() {}
	explicit GF61(const Z61 & s0, const Z61 & s1) : _s0(s0), _s1(s1) {}
	explicit GF61(const uint64_t n0, const uint64_t n1) : _s0(n0), _s1(n1) {}

	const Z61 & s0() const { return _s0; }
	const Z61 & s1() const { return _s1; }

	bool operator==(const GF61 & rhs) const { return ((_s0 == rhs._s0) && (_s1 == rhs._s1)); }

	GF61 conj() const { return GF61(_s0, -_s1); }

	GF61 & operator+=(const GF61 & rhs) { _s0 += rhs._s0; _s1 += rhs._s1; return *this; }
	GF61 & operator-=(const GF61 & rhs) { _s0 -= rhs._s0; _s1 -= rhs._s1; return *this; }
	GF61 & operator*=(const Z61 & rhs) { _s0 *= rhs; _s1 *= rhs; return *this; }
	GF61 operator*=(const GF61 & rhs) { const Z61 t = _s1 * rhs._s0 + _s0 * rhs._s1; _s0 = _s0 * rhs._s0 - _s1 * rhs._s1; _s1 = t; return *this; }

	GF61 operator+(const GF61 & rhs) const { GF61 r = *this; r += rhs; return r; }
	GF61 operator-(const GF61 & rhs) const { GF61 r = *this; r -= rhs; return r; }
	GF61 operator*(const Z61 & rhs) const { GF61 r = *this; r *= rhs; return r; }
	GF61 operator*(const GF61 & rhs) const { GF61 r = *this; r *= rhs; return r; }

	GF61 pow(const __uint128_t e) const
	{
		if (e == 0) return GF61(1u, 0u);
		GF61 r = GF61(1u, 0u), y = *this;
		for (__uint128_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r *= y; y *= y; }
		return r * y;
	}

	GF61 invert() const { __uint128_t t = (uint64_t(1) << 61) - 1; t = t * t - 2; return GF61(pow(t)); }

	static const GF61 root_nth(const size_t n) { return GF61(Z61(_h_0), Z61(_h_1)).pow(_h_order / n); }
};

class GF61v3
{
private:
	GF61 _z[3];

public:
	GF61v3() {}
	explicit GF61v3(const GF61 & z0, const GF61 & z1, const GF61 & z2) { _z[0] = z0; _z[1] = z1; _z[2] = z2; }

	const GF61 & operator[](const size_t i) const { return _z[i]; }

	GF61v3 conj() const { return GF61v3(_z[0].conj(), _z[1].conj(), _z[2].conj()); }

	GF61v3 & operator+=(const GF61v3 & rhs) { for (size_t i = 0; i < 3; ++i) _z[i] += rhs._z[i]; return *this; }
	GF61v3 & operator-=(const GF61v3 & rhs) { for (size_t i = 0; i < 3; ++i) _z[i] -= rhs._z[i]; return *this; }
	GF61v3 & operator*=(const Z61 & rhs) { for (size_t i = 0; i < 3; ++i) _z[i] *= rhs; return *this; }
	GF61v3 & operator*=(const GF61 & rhs) { for (size_t i = 0; i < 3; ++i) _z[i] *= rhs; return *this; }
	GF61v3 & operator*=(const GF61v3 & rhs) { for (size_t i = 0; i < 3; ++i) _z[i] *= rhs._z[i]; return *this; }

	GF61v3 operator+(const GF61v3 & rhs) const { GF61v3 r = *this; r += rhs; return r; }
	GF61v3 operator-(const GF61v3 & rhs) const { GF61v3 r = *this; r -= rhs; return r; }
	GF61v3 operator*(const Z61 & rhs) const { GF61v3 r = *this; r *= rhs; return r; }
	GF61v3 operator*(const GF61 & rhs) const { GF61v3 r = *this; r *= rhs; return r; }
	GF61v3 operator*(const GF61v3 & rhs) { GF61v3 r = *this; r *= rhs; return r; }

	void weight(const GF61 & w)
	{
		_z[1] *= w; _z[2] *= w * w;
	}

	void forward3()
	{
		const GF61 u0 = _z[0], u1 = _z[1], u2 = _z[2];
		const GF61 t = (u1 - u2) * J;
		_z[0] = u0 + u1 + u2;
		_z[1] = u0 - u2 + t;
		_z[2] = u0 - u1 - t;
	}

	void backward3()
	{
		const GF61 u0 = _z[0], u1 = _z[1], u2 = _z[2];
		const GF61 t = (u1 - u2) * J;
		_z[0] = (u0 + u1 + u2) * inv3;
		_z[1] = (u0 - u1 - t) * inv3;
		_z[2] = (u0 - u2 + t) * inv3;
	}
};

// bit-reversal permutation of index i for a sequence of n items
static constexpr size_t bitrev(const size_t i, const size_t n)
{
	size_t r = 0;
	for (size_t k = n, j = i; k > 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
	return r;
}

// Bit-reversed permutation
static void scramble(GF61v3 * const z, const size_t n)
{
	for (size_t i = 0, j = 0; i < n - 1; ++i)
	{
		if (i < j) std::swap(z[i], z[j]);
		size_t k = n / 2; while (k <= j) { j -= k; k /= 2; }
		j += k;
	}
}

static void forward2(GF61v3 * const z, const size_t m, const size_t s)
{
	const GF61 r2m = GF61::root_nth(2 * m);

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF61v3 u0 = z[k + 0 * m], u1 = z[k + 1 * m];
			z[k + 0 * m] = u0 + u1; z[k + 1 * m] = (u0 - u1) * r2m.pow(i);
		}
	}
}

static void backward2(GF61v3 * const z, const size_t m, const size_t s)
{
	const GF61 r2mi = GF61::root_nth(2 * m).conj();

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF61v3 u0 = z[k + 0 * m], u1 = z[k + 1 * m] * r2mi.pow(i);
			z[k + 0 * m] = (u0 + u1) * inv2; z[k + 1 * m] = (u0 - u1) * inv2;
		}
	}
}

// Q = P^2 mod x^n - 1
static void square_slow(uint64_t * const Q, const uint64_t * const P, const size_t n)
{
	for (size_t j = 0; j < n; ++j)
	{
		uint64_t l = 0; for (size_t i = 0; i <= j; ++i) l += P[i] * P[j - i];
		uint64_t h = 0; for (size_t i = j + 1; i < n; ++i) h += P[i] * P[j - i + n];
		Q[j] = l + h;
	}
}

inline void check(const uint64_t * const P, const uint64_t * const Q, const size_t n)
{
	bool error = false;
	uint64_t a_max = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const uint64_t a = P[i], b = Q[i];
		if (a != b) { error = true; std::cout << i << ": " << a << " != " << b << std::endl; }
		a_max = std::max(a_max, a);
	}
	std::cout << (error ? "N" : "") << "OK, max = " << a_max << std::endl;
}

int main()
{
	std::srand((unsigned int)(std::time(nullptr)));

	const size_t n = 3 << 11;
	uint64_t * const P = new uint64_t[n];
	uint64_t * const Q = new uint64_t[n];
	uint64_t * const R = new uint64_t[n];
	for (size_t i = 0; i < n; ++i) P[i] = uint64_t(std::rand()) % 1000u;

	square_slow(Q, P, n);

	GF61v3 * const z = new GF61v3[n / 3];
	for (size_t i = 0; i < n / 3; ++i) z[i] = GF61v3(GF61(P[3 * i + 0], 0), GF61(P[3 * i + 1], 0), GF61(P[3 * i + 2], 0));
	for (size_t m = n / 6, s = 1; m >= 1; m /= 2, s *= 2) forward2(z, m, s);

	const GF61 w_n = GF61::root_nth(n);

	for (size_t i = 0; i < n / 3; ++i)
	{
		const GF61 w_i = w_n.pow(bitrev(i, n / 3));
		z[i].weight(w_i);
		z[i].forward3();
		z[i] *= z[i];
		z[i].backward3();
		z[i].weight(w_i.invert());
	}

	for (size_t m = 1, s = n / 6; s >= 1; m *= 2, s /= 2) backward2(z, m, s);
	for (size_t i = 0; i < n / 3; ++i) for (size_t l = 0; l < 3; ++l) R[3 * i + l] = z[i][l].s0().get();

	check(Q, R, n);

	GF61v3 * const z2 = new GF61v3[n / 6];	// half length

	for (size_t i = 0; i < n / 6; ++i) z2[i] = GF61v3(GF61(P[6 * i + 0], P[6 * i + 3]), GF61(P[6 * i + 1], P[6 * i + 4]), GF61(P[6 * i + 2], P[6 * i + 5]));
	for (size_t m = n / 12, s = 1; m >= 1; m /= 2, s *= 2) forward2(z2, m, s);
	scramble(z2, n / 6);

	const GF61 r = GF61::root_nth(2 * n / 6), ri = r.conj();

	for (size_t i = 0; i < n / 6; ++i)
	{
		const size_t mi = (i == 0) ? 0 : n / 6 - i;
		const GF61v3 z2i = z2[i], z2mi = z2[mi].conj();
		const GF61v3 f = (z2i + z2mi) * inv2, g = (z2i - z2mi) * inv2 * r.pow(i);
		z[i] = f - g * GF61(0, 1);
		if (i == 0) z[n / 6] = f + g * GF61(0, 1);
	}

	for (size_t i = 0; i <= n / 6; ++i)
	{
		const GF61 w_i = w_n.pow(i);
		z[i].weight(w_i);
		z[i].forward3();
		z[i] *= z[i];
		z[i].backward3();
		z[i].weight(w_i.invert());
	}

	for (size_t i = 0; i < n / 6; ++i)
	{
		const GF61v3 zi = z[i], zmi = z[n / 6 - i].conj();
		z2[i] = ((zi + zmi) + (zi - zmi) * ri.pow(i) * GF61(0, 1)) * inv2;
	}

	scramble(z2, n / 6);
	for (size_t m = 1, s = n / 12; s >= 1; m *= 2, s /= 2) backward2(z2, m, s);
	for (size_t i = 0; i < n / 6; ++i) for (size_t l = 0; l < 3; ++l) { R[6 * i + 0 + l] = z2[i][l].s0().get(); R[6 * i + 3 + l] = z2[i][l].s1().get(); }

	check(Q, R, n);

	delete[] z2;
	delete[] z;

	delete[] P;
	delete[] Q;
	delete[] R;

	return EXIT_SUCCESS;
}
