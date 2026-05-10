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
	uint64_t _n;

	static const uint64_t _p = 3 * (uint64_t(1) << 41) + 1;
	static const uint64_t _primroot = 5;

public:
	GF() {}
	GF(const uint64_t n) : _n(n) {}

	uint64_t get() const { return _n; }

	GF & operator+=(const GF & rhs) { const uint64_t c = (_n >= _p - rhs._n) ? _p : 0; _n += rhs._n; _n -= c; return *this; }
	GF & operator-=(const GF & rhs) { const uint64_t c = (_n < rhs._n) ? _p : 0; _n -= rhs._n; _n += c; return *this; }
	GF & operator*=(const GF & rhs) { *this = uint64_t((_n * __uint128_t(rhs._n)) % _p); return *this; }

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

static const GF J = GF::root_nth(3), inv2 = GF(2).invert(), inv3 = GF(3).invert();

static void forward2(GF * const x, const size_t m, const size_t s)
{
	const GF r2m = GF::root_nth(2 * m);

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF u0 = x[k + 0 * m], u1 = x[k + 1 * m];
			x[k + 0 * m] = u0 + u1; x[k + 1 * m] = (u0 - u1) * r2m.pow(i);
		}
	}
}

static void backward2(GF * const x, const size_t m, const size_t s)
{
	const GF r2mi = GF::root_nth(2 * m).invert();

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF u0 = x[k + 0 * m], u1 = x[k + 1 * m] * r2mi.pow(i);
			x[k + 0 * m] = (u0 + u1) * inv2; x[k + 1 * m] = (u0 - u1) * inv2;
		}
	}
}

static void forward3(GF * const x, const size_t m)
{
	for (size_t i = 0; i < m; ++i)
	{
		const GF u0 = x[3 * i + 0], u1 = x[3 * i + 1], u2 = x[3 * i + 2];
		const GF t = J * (u1 - u2);
		x[3 * i + 0] = u0 + u1 + u2;
		x[3 * i + 1] = u0 - u2 + t;
		x[3 * i + 2] = u0 - u1 - t;
	}
}

static void backward3(GF * const x, const size_t m)
{
	for (size_t i = 0; i < m; ++i)
	{
		const GF u0 = x[3 * i + 0], u1 = x[3 * i + 1], u2 = x[3 * i + 2];
		const GF t = J * (u1 - u2);
		x[3 * i + 0] = (u0 + u1 + u2) * inv3;
		x[3 * i + 1] = (u0 - u1 - t) * inv3;
		x[3 * i + 2] = (u0 - u2 + t) * inv3;
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

	const size_t n = 3 << 12;
	uint64_t P[n], Q[n];
	 for (size_t i = 0; i < n; ++i) P[i] = uint64_t(std::rand()) % 1000u;

	square_slow(Q, P, n);

	GF R[n];
	for (size_t i = 0; i < n; ++i) R[i] = GF(P[i]);
	for (size_t m = n / 2, s = 1; m >= 3; m /= 2, s *= 2) forward2(R, m, s);
	forward3(R, n / 3);
	for (size_t k = 0; k < n; ++k) R[k] *= R[k];
	backward3(R, n / 3);
	for (size_t m = 3, s = n / 6; s >= 1; m *= 2, s /= 2) backward2(R, m, s);
	for (size_t i = 0; i < n; ++i) P[i] = R[i].get();

	check(P, Q, n);

	return EXIT_SUCCESS;
}
