/*
Copyright 2025, Yves Gallot

mersenne2.cpp is free source code. You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>
#include <ctime>

// a finite field
class GF
{
private:
	uint64_t n;

	static const uint64_t _p = 3 * (uint64_t(1) << 41) + 1;
	static const uint64_t _primroot = 5;

public:
	GF() {}
	GF(const uint64_t l) : n(l) {}

	uint64_t get() const { return n; }

	GF & operator+=(const GF & rhs) { const uint64_t c = (n >= _p - rhs.n) ? _p : 0; n += rhs.n; n -= c; return *this; }
	GF & operator-=(const GF & rhs) { const uint64_t c = (n < rhs.n) ? _p : 0; n -= rhs.n; n += c; return *this; }
	GF & operator*=(const GF & rhs) { *this = uint64_t((n * __uint128_t(rhs.n)) % _p); return *this; }

	GF operator+(const GF & rhs) const { GF r = *this; r += rhs; return r; }
	GF operator-(const GF & rhs) const { GF r = *this; r -= rhs; return r; }
	GF operator*(const GF & rhs) const { GF r = *this; r *= rhs; return r; }

	GF pow(const size_t e) const
	{
		if (e == 0) return GF(1u);

		GF r = GF(1u), y = *this;
		for (size_t i = e; i != 1; i /= 2)
		{
			if (i % 2 != 0) r *= y;
			y *= y;
		}
		r *= y;

		return r;
	}

	GF invert() const { return pow(_p - 2); }

	static const GF root_nth(const size_t n) { return GF(_primroot).pow((_p - 1) / n); }
};

// bit-reversal permutation of index i for a sequence of n items
static constexpr size_t bitrev(const size_t i, const size_t n)
{
	size_t r = 0;
	for (size_t k = n, j = i; k > 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
	return r;
}

static void forward2(GF * const z, const size_t m, const size_t s)
{
	const GF ws = GF::root_nth(2 * s);

	for (size_t j = 0; j < s; ++j)
	{
		const GF wj = ws.pow(bitrev(j, s));
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF u0 = z[k + 0 * m], u1 = z[k + 1 * m] * wj;
			z[k + 0 * m] = u0 + u1; z[k + 1 * m] = u0 - u1;
		}
	}
}

static void backward2(GF * const z, const size_t m, const size_t s)
{
	const GF wsi = GF::root_nth(2 * s).invert(), inv2 = GF(2).invert();

	for (size_t j = 0; j < s; ++j)
	{
		const GF wji = wsi.pow(bitrev(j, s));

		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF u0 = z[k + 0 * m], u1 = z[k + 1 * m];
			z[k + 0 * m] = (u0 + u1) * inv2;
			z[k + 1 * m] = (u0 - u1) * wji * inv2;
		}
	}
}

static void forward3(GF * const z, const size_t m)
{
	static const GF J = GF::root_nth(3), J2 = J * J;
	const GF rm = GF::root_nth(3 * m);

	for (size_t i = 0; i < m; ++i)
	{
		const GF u0 = z[i + 0 * m], u1 = z[i + 1 * m], u2 = z[i + 2 * m];
		z[i + 0 * m] = (u0 +       u1 +      u2);
		z[i + 1 * m] = (u0 +   J * u1 + J2 * u2) * rm.pow(1 * i);
		z[i + 2 * m] = (u0 +  J2 * u1 +  J * u2) * rm.pow(2 * i);
	}
}

static void backward3(GF * const z, const size_t m)
{
	static const GF J = GF::root_nth(3), J2 = J * J, inv3 = GF(3).invert();
	const GF rmi = GF::root_nth(3 * m).invert();

	for (size_t i = 0; i < m; ++i)
	{
		const GF u0 = z[i + 0 * m], u1 = z[i + 1 * m] * rmi.pow(1 * i), u2 = z[i + 2 * m] * rmi.pow(2 * i);
		z[i + 0 * m] = (u0 +      u1 +      u2) * inv3;
		z[i + 1 * m] = (u0 + J2 * u1 +  J * u2) * inv3;
		z[i + 2 * m] = (u0 +  J * u1 + J2 * u2) * inv3;
	}
}

static void sqr(GF * const z, const size_t n)
{
	for (size_t k = 0; k < n; ++k) z[k] *= z[k];
}

static void square2(GF * const z, const size_t n)
{
	for (size_t m = n / 2, s = 1; m >= 1; m /= 2, s *= 2) forward2(z, m, s);
	sqr(z, n);
	for (size_t m = 1, s = n / 2; s >= 1; m *= 2, s /= 2) backward2(z, m, s);
}

// Q = P^2 mod x^n - 1
static void square_slow(GF * const Q, const GF * const P, const size_t n)
{
	for (size_t j = 0; j < n; ++j)
	{
		GF l = 0; for (size_t i = 0; i <= j; ++i) l += P[i] * P[j - i];
		GF h = 0; for (size_t i = j + 1; i < n; ++i) h += P[i] * P[j - i + n];
		Q[j] = l + h;
	}
}

inline void check(const GF * const P, const GF * const Q, const size_t n)
{
	bool error = false;
	uint64_t a_max = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const uint64_t a = P[i].get(), b = Q[i].get();
		if (a != b) { error = true; std::cout << i << ": " << a << " != " << b << std::endl; }
		a_max = std::max(a_max, a);
	}
	std::cout << (error ? "N" : "") << "OK, max = " << a_max << std::endl;
}

int main()
{
	std::srand((unsigned int)(std::time(nullptr)));

	const size_t m = 1 << 12, n = 3 * m;
	GF P[n]; for (size_t i = 0; i < n; ++i) P[i] = GF(uint64_t(std::rand()) % 1000u);

	GF Q[n]; square_slow(Q, P, n);

	forward3(P, m);
	square2(&P[0 * m], m); square2(&P[1 * m], m); square2(&P[2 * m], m);
	backward3(P, m);

	check(P, Q, n);

	return EXIT_SUCCESS;
}
