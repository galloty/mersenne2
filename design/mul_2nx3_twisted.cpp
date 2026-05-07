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

static void square3(GF * const P)
{
	static const GF J = GF::root_nth(3), J2 = J * J, inv3 = GF(3).invert();

	const GF u0 = P[0], u1 = P[1], u2 = P[2];
	GF v0 = (u0 +       u1 +      u2);
	GF v1 = (u0 +   J * u1 + J2 * u2);
	GF v2 = (u0 +  J2 * u1 +  J * u2);
	v0 *= v0; v1 *= v1; v2 *= v2;
	P[0] = (v0 +      v1 +      v2) * inv3;
	P[1] = (v0 + J2 * v1 +  J * v2) * inv3;
	P[2] = (v0 +  J * v1 + J2 * v2) * inv3;

	// const GF c = P[0], b = P[1], a = P[2];
	// P[0] = (a * b) * 2 + c * c;
	// P[1] = (b * c) * 2 + a * a;
	// P[2] = (c * a) * 2 + b * b;
}

static void forward2(GF * const z, const size_t m, const size_t s)
{
	const GF wm = GF::root_nth(2 * m);

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF u0 = z[k + 0 * m], u1 = z[k + 1 * m];
			z[k + 0 * m] = u0 + u1; z[k + 1 * m] = (u0 - u1) * wm.pow(i);;
		}
	}
}

static void backward2(GF * const z, const size_t m, const size_t s)
{
	const GF wmi = GF::root_nth(2 * m).invert(), inv2 = GF(2).invert();

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const GF u0 = z[k + 0 * m], u1 = z[k + 1 * m] * wmi.pow(i);
			z[k + 0 * m] = (u0 + u1) * inv2; z[k + 1 * m] = (u0 - u1) * inv2;
		}
	}
}

static void sqr(GF * const z, const size_t m)
{
	for (size_t i = 0; i < m; ++i) square3(&z[3 * i]);
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

	const size_t n = 3 << 12;
	GF P[n]; for (size_t i = 0; i < n; ++i) P[i] = GF(uint64_t(std::rand()) % 1000u);

	GF Q[n]; square_slow(Q, P, n);

	for (size_t m = n / 2, s = 1; m >= 3; m /= 2, s *= 2) forward2(P, m, s);
	sqr(P, n / 3);
	for (size_t m = 3, s = n / 6; s >= 1; m *= 2, s /= 2) backward2(P, m, s);

	check(P, Q, n);

	return EXIT_SUCCESS;
}
