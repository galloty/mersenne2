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

static const GF J = GF::root_nth(3), inv3 = GF(3).invert();

// a quadratic polynomial modulo x^3 - 1 over the finite field
class Quad
{
private:
	GF _x[3];

public:
	Quad() {}
	Quad(const GF x[3]) { for (size_t i = 0; i < 3; ++i) _x[i] = x[i]; }

	void get(GF x[3]) const { for (size_t i = 0; i < 3; ++i) x[i] = _x[i]; }

	Quad & operator+=(const Quad & rhs) { for (size_t i = 0; i < 3; ++i) _x[i] += rhs._x[i]; return *this; }
	Quad & operator-=(const Quad & rhs) { for (size_t i = 0; i < 3; ++i) _x[i] -= rhs._x[i]; return *this; }
	Quad & operator*=(const GF & rhs) { for (size_t i = 0; i < 3; ++i) _x[i] *= rhs; return *this; }

	Quad operator+(const Quad & rhs) const { Quad r = *this; r += rhs; return r; }
	Quad operator-(const Quad & rhs) const { Quad r = *this; r -= rhs; return r; }
	Quad operator*(const GF & rhs) const { Quad r = *this; r *= rhs; return r; }

	Quad dot(const Quad & rhs) const { Quad r = *this; for (size_t i = 0; i < 3; ++i) r._x[i] *= rhs._x[i]; return r; }

	void forward3()
	{
		const GF u0 = _x[0], u1 = _x[1], u2 = _x[2];
		const GF t = J * (u1 - u2);
		_x[0] = u0 + u1 + u2;
		_x[1] = u0 - u2 + t;
		_x[2] = u0 - u1 - t;
	}

	void backward3()
	{
		const GF u0 = _x[0], u1 = _x[1], u2 = _x[2];
		const GF t = J * (u1 - u2);
		_x[0] = u0 + u1 + u2;
		_x[1] = u0 - u1 - t;
		_x[2] = u0 - u2 + t;
		*this *= inv3;
	}
};

inline Quad pow_w(const GF & w, const size_t i)
{
	GF wi[3]; wi[0] = w.pow(3 * i + 0); wi[1] = wi[0] * w; wi[2] = wi[1] * w;
	return Quad(wi);
}

static void forward2(Quad * const q, const size_t m, const size_t s)
{
	const GF wm = GF::root_nth(6 * m);

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const Quad u0 = q[k + 0 * m], u1 = q[k + 1 * m];
			q[k + 0 * m] = u0 + u1; q[k + 1 * m] = (u0 - u1).dot(pow_w(wm, i));
		}
	}
}

static void backward2(Quad * const q, const size_t m, const size_t s)
{
	const GF wmi = GF::root_nth(6 * m).invert(), inv2 = GF(2).invert();

	for (size_t j = 0; j < s; ++j)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const size_t k = 2 * m * j + i;
			const Quad u0 = q[k + 0 * m], u1 = q[k + 1 * m].dot(pow_w(wmi, i));
			q[k + 0 * m] = (u0 + u1) * inv2; q[k + 1 * m] = (u0 - u1) * inv2;
		}
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

	GF R[n]; square_slow(R, P, n);

	Quad Q[n / 3];
	for (size_t i = 0; i < n / 3; ++i) Q[i] = Quad(&P[3 * i]);

	for (size_t m = n / 6, s = 1; m >= 1; m /= 2, s *= 2) forward2(Q, m, s);
	for (size_t i = 0; i < n / 3; ++i)
	{
		Quad & Qi = Q[i];
		Qi.forward3();
		Qi = Qi.dot(Qi);
		Qi.backward3();
	}
	for (size_t m = 1, s = n / 6; s >= 1; m *= 2, s /= 2) backward2(Q, m, s);

	for (size_t i = 0; i < n / 3; ++i) Q[i].get(&P[3 * i]);

	check(P, R, n);

	return EXIT_SUCCESS;
}
