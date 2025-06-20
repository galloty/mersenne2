/*
Copyright 2025, Yves Gallot

From: Nick Craig-Wood, IOCCC 2012 Entry, https://github.com/ncw/ioccc2012
NTT algorithm is Gentleman-Sande and its inverse is Cooley-Tukey recursion.
Transform length is 2^m or 5*2^m. Butterflies are radix-4 and radix-5.

mersenne_NCW45.cpp is free source code. You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <iostream>
#include <cstdint>
#include <cmath>

class Zp
{
private:
	uint64_t n;

public:
	static const uint64_t p = (((1ull << 32) - 1) << 32) + 1;	// 2^64 - 2^32 + 1

private:
	Zp reduce(const __uint128_t t) const
	{
		const uint64_t hi = uint64_t(t >> 64), lo = uint64_t(t);

		// hih * 2^96 + hil * 2^64 + lo = lo + hil * 2^32 - (hih + hil)
		Zp r = Zp((lo >= p) ? lo - p : lo);
		r += Zp(hi << 32);				// lhs * rhs < p^2 => hi * 2^32 < p^2 / 2^32 < p.
		r -= Zp((hi >> 32) + uint32_t(hi));
		return r;
	}

public:
	Zp() {}
	Zp(const uint64_t l) : n(l) {}

	uint64_t get() const { return n; }
	void set(const uint64_t l) { n = l; }

	Zp & operator+=(const Zp & rhs) { const uint64_t c = (n >= p - rhs.n) ? p : 0; n += rhs.n; n -= c; return *this; }
	Zp & operator-=(const Zp & rhs) { const uint64_t c = (n < rhs.n) ? p : 0; n -= rhs.n; n += c; return *this; }
	Zp & operator*=(const Zp & rhs) { *this = reduce(n * __uint128_t(rhs.n)); return *this; }

	Zp operator+(const Zp & rhs) const { Zp r = *this; r += rhs; return r; }
	Zp operator-(const Zp & rhs) const { Zp r = *this; r -= rhs; return r; }
	Zp operator*(const Zp & rhs) const { return reduce(n * __uint128_t(rhs.n)); }

	Zp muli() const { return reduce(__uint128_t(n) << 48); }	// sqrt(-1) = 2^48 (mod p)

	Zp pow(const uint64_t e) const
	{
		if (e == 0) return Zp(1u);

		Zp r = Zp(1u), y = *this;
		for (uint64_t i = e; i != 1; i /= 2)
		{
			if (i % 2 != 0) r *= y;
			y *= y;
		}
		r *= y;

		return r;
	}

	Zp invert() const { return Zp(pow(p - 2)); }
	static const Zp primroot(const size_t n) { return Zp(7u).pow((p - 1) / n); }
};

class mersenne
{
private:
	const size_t _n, _n5;
	Zp * const _x;
	Zp * const _w;
	Zp * const _invw;
	Zp * const _w5;
	Zp * const _invw5;
	Zp * const _digit_weight;
	Zp * const _digit_invweight;
	int * const _digit_width;
	bool _even_exponent;

private:
	static constexpr size_t transformsize(const uint32_t exponent)
	{
		// Make sure the transform is long enough so that each 'digit' can't overflow after the convolution.
		uint32_t w = 0, log2_n = 0, log2_n5 = 0;
		do
		{
			++log2_n;
			// digit-width is w or w + 1
			w = exponent >> log2_n;
		// The condition is n * (2^{w + 1} - 1)^2 < 2^64 - 2^32 + 1.
		// If (w + 1) * 2 + log2(n) = 63 then n * (2^{w + 1} - 1)^2 < n * (2^{w + 1})^2 = 2^63 < 2^64 - 2^32 + 1.
		} while ((w + 1) * 2 + log2_n >= 64);

		do
		{
			++log2_n5;
			w = exponent / (5u << log2_n5);
		// log2(5) ~ 2.3219 < 2.4
		} while ((w + 1) * 2 + (log2_n5 + 2.4) >= 64);

		return std::min(size_t(1) << log2_n, size_t(5) << log2_n5);
	}

	// Add a carry onto the number and return the carry of the first digit_width bits
	static constexpr uint64_t digit_adc(const uint64_t lhs, const int digit_width, uint64_t & carry)
	{
		const uint64_t s = lhs + carry;
		const uint64_t c =  s < lhs;
		carry = (s >> digit_width) + (c << (64 - digit_width));
		return s & ((uint64_t(1) << digit_width) - 1);
	}

	void mul(const Zp * const y) const
	{
		Zp * const x = _x;
		for (size_t k = 0, n = _n; k < n; ++k) x[k] *= y[k];
	}

	void ntt() const
	{
		const size_t n = _n, n5 = _n5;
		Zp * const x = _x;
		const Zp * const w = _w;
		const Zp * const w5 = _w5;

		if (n % 5 == 0)
		{
			static const Zp r1 = Zp::primroot(5u), r2 = r1 * r1, r3 = r2 * r1, r4 = r2 * r2;

			// Radix-5
			for (size_t k = 0; k < n5; ++k)
			{
				const Zp t0 = x[k + 0 * n5], t1 = x[k + 1 * n5], t2 = x[k + 2 * n5], t3 = x[k + 3 * n5], t4 = x[k + 4 * n5];
				x[k + 0 * n5] =  t0 +      t1 +      t2 +      t3 +      t4;
				x[k + 1 * n5] = (t0 + r1 * t1 + r2 * t2 + r3 * t3 + r4 * t4) * w5[4 * k + 0];
				x[k + 2 * n5] = (t0 + r2 * t1 + r4 * t2 + r1 * t3 + r3 * t4) * w5[4 * k + 1];
				x[k + 3 * n5] = (t0 + r3 * t1 + r1 * t2 + r4 * t3 + r2 * t4) * w5[4 * k + 2];
				x[k + 4 * n5] = (t0 + r4 * t1 + r3 * t2 + r2 * t3 + r1 * t4) * w5[4 * k + 3];
			}
		}

		// Radix-4
		for (size_t m = n5 / 4; m >= 1; m /= 4)
		{
			const Zp * const wm = &w[3 * 2 * m];

			for (size_t k = 0; k < n / 4; ++k)
			{
				const size_t j = k & (m - 1), i = 4 * (k - j) + j;

				const Zp w2 = wm[3 * j + 0], w1 = wm[3 * j + 1], w12 = wm[3 * j + 2];

				const Zp u0 = x[i + 0 * m], u1 = x[i + 1 * m], u2 = x[i + 2 * m], u3 = x[i + 3 * m];
				const Zp v0 = u0 + u2, v1 = u1 + u3, v2 = u0 - u2, v3 = Zp(u1 - u3).muli();
				x[i + 0 * m] = v0 + v1;
				x[i + 1 * m] = (v0 - v1) * w1;
				x[i + 2 * m] = (v2 + v3) * w2;
				x[i + 3 * m] = (v2 - v3) * w12;
			}
		}
	}

	void nttinv() const
	{
		const size_t n = _n, n5 = _n5;
		Zp * const x = _x;
		const Zp * const invw = _invw;
		const Zp * const invw5 = _invw5;

		// Radix-4
		for (size_t m = _even_exponent ? 1 : 2; m <= n5 / 4; m *= 4)
		{
			const Zp * const invwm = &invw[3 * 2 * m];

			for (size_t k = 0; k < n / 4; ++k)
			{
				const size_t j = k & (m - 1), i = 4 * (k - j) + j;

				const Zp iw2 = invwm[3 * j + 0], iw1 = invwm[3 * j + 1], iw12 = invwm[3 * j + 2];

				const Zp u0 = x[i + 0 * m], u1 = x[i + 1 * m] * iw1, u2 = x[i + 2 * m] * iw2, u3 = x[i + 3 * m] * iw12;
				const Zp v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp(u3 - u2).muli();
				x[i + 0 * m] = v0 + v2;
				x[i + 1 * m] = v1 + v3;
				x[i + 2 * m] = v0 - v2;
				x[i + 3 * m] = v1 - v3;
			}
		}

		if (n % 5 == 0)
		{
			static const Zp r1 = Zp::primroot(5u), r2 = r1 * r1, r3 = r2 * r1, r4 = r2 * r2;

			// Radix-5
			for (size_t k = 0; k < n5; ++k)
			{
				const Zp t0 = x[k + 0 * n5], t1 = x[k + 1 * n5] * invw5[4 * k + 0], t2 = x[k + 2 * n5] * invw5[4 * k + 1];
				const Zp t3 = x[k + 3 * n5] * invw5[4 * k + 2], t4 = x[k + 4 * n5] * invw5[4 * k + 3];
				x[k + 0 * n5] = t0 + t1 + t2 + t3 + t4;
				x[k + 1 * n5] = t0 + r4 * t1 + r3 * t2 + r2 * t3 + r1 * t4;
				x[k + 2 * n5] = t0 + r3 * t1 + r1 * t2 + r4 * t3 + r2 * t4;
				x[k + 3 * n5] = t0 + r2 * t1 + r4 * t2 + r1 * t3 + r3 * t4;
				x[k + 4 * n5] = t0 + r1 * t1 + r2 * t2 + r3 * t3 + r4 * t4;
			}
		}
	}

	void square1() const
	{
		const size_t n = _n;
		Zp * const x = _x;

		for (size_t k = 0; k < n; ++k) x[k] *= x[k];
	}

	void square2() const
	{
		const size_t n = _n;
		Zp * const x = _x;

		for (size_t k = 0; k < n; k += 2)
		{
			// Radix-2 NTT, square, Radix-2 inverse NTT
			const Zp u0 = x[k + 0], u1 = x[k + 1];
			const Zp v0 = u0 + u1, v1 = u0 - u1;
			const Zp s0 = v0 * v0, s1 = v1 * v1;
			x[k + 0] = s0 + s1; x[k + 1] = s0 - s1;
		}
	}

	void carry() const
	{
		const size_t n = _n;
		Zp * const x = _x;
		const int * const digit_width = _digit_width;

		uint64_t c = 0;
		for (size_t k = 0; k < n; ++k) x[k].set(digit_adc(x[k].get(), digit_width[k], c));

		while (c != 0)
		{
			for (size_t k = 0; k < n; ++k)
			{
				x[k].set(digit_adc(x[k].get(), digit_width[k], c));
				if (c == 0) break;
			}
		}
	}

public:
	mersenne(const uint32_t q) : _n(transformsize(q)), _n5((_n % 5 == 0) ? _n / 5 : _n),
		_x(new Zp[_n]), _w(new Zp[3 * _n5]), _invw(new Zp[3 * _n5]),
		_w5((_n % 5 == 0) ? new Zp[4 * _n5] : nullptr), _invw5((_n % 5 == 0) ? new Zp[4 * _n5] : nullptr),
		_digit_weight(new Zp[_n]), _digit_invweight(new Zp[_n]), _digit_width(new int[_n])
	{
		const size_t n = _n, n5 = _n5;

		size_t m = n5; for (; m > 1; m /= 4);
		_even_exponent = (m == 1);

		// twiddle factors

		// Radix-4
		Zp * const w = _w;
		Zp * const invw = _invw;

		for (size_t m = n5 / 2, s = 1; m >= 1; m /= 2, s *= 2)
		{
			const Zp r = Zp::primroot(2 * m), invr = r.invert();
			Zp w_j = Zp(1u), invw_j = Zp(1u);
			for (size_t j = 0; j < m; j++)
			{
				w[3 * (m + j) + 0] = w_j;
				const Zp w_j2 = w_j * w_j;
				w[3 * (m + j) + 1] = w_j2;
				w[3 * (m + j) + 2] = w_j2 * w_j;

				invw[3 * (m + j) + 0] = invw_j;
				const Zp invw_j2 = invw_j * invw_j;
				invw[3 * (m + j) + 1] = invw_j2;
				invw[3 * (m + j) + 2] = invw_j2 * invw_j;

				w_j *= r; invw_j *= invr;
			}
		}

		if (n % 5 == 0)
		{
			// Radix-5
			Zp * const w5 = _w5;
			Zp * const invw5 = _invw5;

			const Zp r = Zp::primroot(n), invr = r.invert();
			Zp w_j = Zp(1u), invw_j = Zp(1u);
			for (size_t j = 0; j < n5; ++j)
			{
				w5[4 * j + 0] = w_j;
				const Zp w_j2 = w_j * w_j;
				w5[4 * j + 1] = w_j2;
				w5[4 * j + 2] = w_j * w_j2;
				w5[4 * j + 3] = w_j2 * w_j2;

				invw5[4 * j + 0] = invw_j;
				const Zp invw_j2 = invw_j * invw_j;
				invw5[4 * j + 1] = invw_j2;
				invw5[4 * j + 2] = invw_j * invw_j2;
				invw5[4 * j + 3] = invw_j2 * invw_j2;

				w_j *= r; invw_j *= invr;
			}
		}

		// weights
		Zp * const digit_weight = _digit_weight;
		Zp * const digit_invweight = _digit_invweight;
		int * const digit_width = _digit_width;

		// n-th root of two
		const Zp nr2 = Zp(554u).pow((Zp::p - 1) / 192 / n);
		const Zp inv_n = Zp(n).invert();

		const uint32_t q_n = q / n;

		digit_weight[0] = Zp(1u);
		digit_invweight[0] = inv_n;

		uint32_t o = 0;
		for (size_t j = 0; j <= n; ++j)
		{
			const uint64_t qj = q * uint64_t(j);
			// ceil(a / b) = floor((a - 1) / b) + 1
			const uint32_t ceil_qj_n = (qj == 0) ? 0 : uint32_t((qj - 1) / n + 1);

			// bit position for digit[i] is ceil(qj / n)
			if (j > 0)
			{
				const uint32_t c = ceil_qj_n - o;
				if ((c != q_n) && (c != q_n + 1)) throw;
				digit_width[j - 1] = c;

				// weight is 2^[ceil(qj / n) - qj/n]
				if (j < n)
				{
					// e = (ceil(qj / n).n - qj) / n
					// qj = k.n => e = 0
					// qj = k.n + r, r > 0 => ((k + 1).n - k.n + r) / n = (n - r) / n
					const uint32_t r = uint32_t(qj % n);
					const Zp nr2r = (r != 0) ? nr2.pow(n - r) : Zp(1u);
					digit_weight[j] = nr2r;
					digit_invweight[j] = nr2r.invert() * inv_n;
				}
			}

			o = ceil_qj_n;
		}
	}

	virtual ~mersenne()
	{
		delete[] _x;
		delete[] _w;
		delete[] _invw;
		if (_n % 5 == 0)
		{
			delete[] _w5;
			delete[] _invw5;
		}
		delete[] _digit_weight;
		delete[] _digit_invweight;
		delete[] _digit_width;
	}

	size_t get_length() const { return _n; }

	void init(const size_t a) const
	{
		const size_t n = _n;
		Zp * const x = _x;

		x[0] = Zp(a);
		for (size_t k = 1; k < n; ++k) x[k] = Zp(0u);
	}

	void square() const
	{
		// weighted convolution
		mul(_digit_weight);
		ntt();
		if (_even_exponent) square1(); else square2();
		nttinv();
		mul(_digit_invweight);

		// carry propagation
		carry();
	}

	void sub(const uint32_t a) const
	{
		const size_t n = _n;
		Zp * const x = _x;
		const int * const digit_width = _digit_width;

		uint32_t c = a;
		while (c != 0)
		{
			for (size_t k = 0; k < n; ++k)
			{
				const uint32_t b = uint32_t(1) << digit_width[k];
				const uint64_t xk = x[k].get();
				if (xk >= c)
				{
					x[k].set(xk - c);
					return;	// done if no carry
				}
				x[k].set(xk - c + b);
				c = 1;
			}
		}
	}

	bool is_zero() const
	{
		const size_t n = _n;
		Zp * const x = _x;

		for (size_t k = 0; k < n; ++k) if (x[k].get() != 0) return false;
		return true;
	}

	bool is_Mp() const
	{
		const size_t n = _n;
		Zp * const x = _x;
		const int * const digit_width = _digit_width;

		for (size_t k = 0; k < n; ++k) if (x[k].get() != (uint64_t(1) << digit_width[k]) - 1) return false;
		return true;
	}
};

int main()
{
	size_t count = 0, count5 = 0;

	// 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, ...
	for (uint32_t p = 3; p <= 1207959503; p += 2)
	{
		bool isprime = true;
		for (uint32_t d = 3; p / d >= d; d += 2) if (p % d == 0) { isprime = false; break; }
		if (!isprime) continue;

		mersenne m(p);

		++count;
		if (m.get_length() % 5 == 0) ++count5;

		m.init(4);
		for (uint32_t i = 0; i < p - 2; ++i)
		{
			m.square();
			m.sub(2);
		}

		if (m.is_zero() || m.is_Mp()) std::cout << p << ": " << m.get_length() << " (radix-5: " << count5 << "/" << count << " = " << 100.0 * count5 / count <<  "%)" << std::endl;
	}

	return EXIT_SUCCESS;
}
