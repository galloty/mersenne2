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

static const GF J = GF::root_nth(3), inv2 = GF(2).invert(), inv3 = GF(3).invert();

static void square3(GF * const P)
{
	static GF J2 = J * J;

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

static void square2twisted(GF * const P, const size_t m)
{
	if (m == 3) { square3(P); return; }		// last stage

	const GF rm = GF::root_nth(m), rmi = rm.invert();
 
	GF * const P0 = &P[0 * m / 2]; GF * const P1 = &P[1 * m / 2];

	for (size_t i = 0; i < m / 2; ++i)
	{
		const GF u0 = P0[i], u1 = P1[i];
		P0[i] = u0 + u1;
		// R1 = x^{m/2} + r, P1 must be twisted such that R1' = x^{m/2} - r
		P1[i] = (u0 - u1) * rm.pow(i);
	}

	square2twisted(P0, m / 2); square2twisted(P1, m / 2);

	for (size_t i = 0; i < m / 2; ++i)
	{
		const GF u0 = P0[i], u1 = P1[i] * rmi.pow(i);	// untwist;
		P0[i] = (u1 + u0) * inv2;
		P1[i] = (u0 - u1) * inv2;
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

	const size_t n = 3 * (1 << 12);
	uint64_t P[n], Q[n];
	for (size_t i = 0; i < n; ++i) P[i] = uint64_t(std::rand()) % 1000u;

	square_slow(Q, P, n);

	GF R[n];
	for (size_t i = 0; i < n; ++i) R[i] = GF(P[i]);
	square2twisted(R, n);
	for (size_t i = 0; i < n; ++i) P[i] = R[i].get();

	check(P, Q, n);

	return EXIT_SUCCESS;
}
