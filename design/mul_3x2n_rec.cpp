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

static void square2(GF * const P, const size_t m, const size_t s, const size_t j)
{
	if (m == 1) { P[0] *= P[0]; return; } 	// P is a scalar

	const GF r = GF::root_nth(2 * s).pow(j), ri = r.invert(), inv2 = GF(2).invert();

	GF * const P0 = &P[0 * m / 2]; GF * const P1 = &P[1 * m / 2];

	for (size_t i = 0; i < m / 2; ++i)
	{
		const GF u0 = P0[i], u1 = r * P1[i];
		P0[i] = u0 + u1; P1[i] = u0 - u1;
	}

	square2(P0, m / 2, s * 2, j + 0 * s); square2(P1, m / 2, s * 2, j + 1 * s);

	for (size_t i = 0; i < m / 2; ++i)
	{
		const GF u0 = P0[i], u1 = P1[i];
		P0[i] = (u1 + u0) * inv2; P1[i] = (u0 - u1) * ri * inv2;
	}
}

static void square3twisted(GF * const P, const size_t m)
{
	static const GF J = GF::root_nth(3), J2 = J * J, inv3 = GF(3).invert();
	const GF rm = GF::root_nth(m), rmi = rm.invert();

	GF * const P0 = &P[0 * m / 3]; GF * const P1 = &P[1 * m / 3]; GF * const P2 = &P[2 * m / 3];

	for (size_t i = 0; i < m / 3; ++i)
	{
		const GF u0 = P0[i], u1 = P1[i], u2 = P2[i];
		P0[i] = (u0 +       u1 +      u2);
		P1[i] = (u0 +   J * u1 + J2 * u2) * rm.pow(1 * i);	// twist
		P2[i] = (u0 +  J2 * u1 +  J * u2) * rm.pow(2 * i);	// twist
	}

	square2(P0, m / 3, 3, 0); square2(P1, m / 3, 3, 0); square2(P2, m / 3, 3, 0);

	for (size_t i = 0; i < m / 3; ++i)
	{
		const GF u0 = P0[i], u1 = P1[i] * rmi.pow(1 * i), u2 = P2[i] * rmi.pow(2 * i);	// untwist
		P0[i] = (u0 +      u1 +      u2) * inv3;
		P1[i] = (u0 + J2 * u1 +  J * u2) * inv3;
		P2[i] = (u0 +  J * u1 + J2 * u2) * inv3;
	}
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

	square3twisted(P, n);

	check(P, Q, n);

	return EXIT_SUCCESS;
}
