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
	uint64_t _n;

	static const uint64_t _p = (((1ull << 32) - 1) << 32) + 1;	// 2^64 - 2^32 + 1
	static const uint32_t _mp64 = uint32_t(-1);					// -p mod (2^64) = 2^32 - 1

private:
	Zp reduce(const __uint128_t t) const
	{
		const uint64_t hi = uint64_t(t >> 64), lo = uint64_t(t);
		const uint32_t hi_lo = uint32_t(hi), hi_hi = uint32_t(hi >> 32);

		// Let X = hi_lo * (2^31 - 1) - hi_hi + lo = hi_hi * 2^96 + hi_lo * 2^64 + lo (mod p)
		// The trick is to add 2^31 - 1 to X (Nick Craig-Wood, ARM-32 assembly code)
		const uint32_t d = _mp64 - hi_hi;
		const uint64_t s = ((uint64_t(hi_lo) << 32) - hi_lo) + d;	// No carry: 0 <= s <= (2^32 - 1)^2 + 2^32 - 1 = 2^64 - 2^32 < p
		uint64_t r = s + lo;	// If carry then r + 2^64 = X + 2^32 - 1. We have r = X (mod p) and r < s < p
		if (r >= s)				// No carry
		{
			// Subtract 2^31 - 1. If the difference is negative then add p (+p = -(-p))
			const uint32_t c = (r < _mp64) ? _mp64 : 0;	// borrow
			r -= _mp64; r -= c;
		}
		return Zp(r);
	}

public:
	Zp() {}
	Zp(const uint64_t n) : _n(n) {}

	static uint64_t p() { return _p; }

	uint64_t get() const { return _n; }
	void set(const uint64_t n) { _n = n; }

	bool operator!=(const Zp & rhs) const { return (_n != rhs._n); }

	Zp & operator+=(const Zp & rhs) { const uint32_t c = (_n >= _p - rhs._n) ? _mp64 : 0; _n += rhs._n; _n += c; return *this; }
	Zp & operator-=(const Zp & rhs) { const uint32_t c = (_n < rhs._n) ? _mp64 : 0; _n -= rhs._n; _n -= c; return *this; }
	Zp & operator*=(const Zp & rhs) { *this = reduce(_n * __uint128_t(rhs._n)); return *this; }

	Zp operator+(const Zp & rhs) const { Zp r = *this; r += rhs; return r; }
	Zp operator-(const Zp & rhs) const { Zp r = *this; r -= rhs; return r; }
	Zp operator*(const Zp & rhs) const { return reduce(_n * __uint128_t(rhs._n)); }

	Zp muli() const { return reduce(__uint128_t(_n) << 48); }	// sqrt(-1) = 2^48 (mod p)

	Zp half() const { return Zp((_n % 2 == 0) ? _n / 2 : ((_n - 1) / 2 + (_p + 1) / 2)); }

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

	Zp invert() const { return Zp(pow(_p - 2)); }
	static const Zp primroot(const size_t n) { return Zp(7u).pow((_p - 1) / n); }
};

class Mersenne
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
	static constexpr uint32_t digit_adc(const uint64_t lhs, const int digit_width, uint64_t & carry)
	{
		const uint64_t s = lhs + carry;
		const uint64_t c = (s < lhs) ? 1 : 0;
		carry = (s >> digit_width) + (c << (64 - digit_width));
		return uint32_t(s) & ((uint32_t(1) << digit_width) - 1);
	}

	void ntt(Zp * const x) const
	{
		const size_t n = _n, n5 = _n5;
		const Zp * const w = _w;
		const Zp * const w5 = _w5;

		if (n % 5 == 0)
		{
			static const Zp K = Zp::primroot(5u), K1 = K + Zp(1u), K2 = K * K;
			static const Zp F2 = (K1 * K2).half(), F3 = K * K1, F4 = K * (K2 + Zp(1u)), F5 = (F3 + F4).half();

			// Radix-5
			for (size_t k = 0; k < n5; ++k)
			{
				const Zp u0 = x[k + 0 * n5], u1 = x[k + 1 * n5], u2 = x[k + 2 * n5], u3 = x[k + 3 * n5], u4 = x[k + 4 * n5];

				// 20 mul, 20 add
				// x[k + 0 * n5] = u0 + u1 + u2 + u3 + u4;
				// x[k + 1 * n5] = (u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4) * w5[4 * k + 0];
				// x[k + 2 * n5] = (u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4) * w5[4 * k + 1];
				// x[k + 3 * n5] = (u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4) * w5[4 * k + 2];
				// x[k + 4 * n5] = (u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4) * w5[4 * k + 3];

				// 8 mul, 20 add
				const Zp v1 = u1 + u4, v4 = u1 - u4, v2 = u2 + u3, v3 = u2 - u3;
				const Zp z0 = v1 + v2;
				const Zp x1 = u0, x2 = (v1 - v2) * F2, x3 = v3 * F3, x4 = v4 * F4, x5 = (v4 - v3) * F5;
				const Zp y1 = x1 + x2, y2 = x1 - x2, y3 = x3 + x5, y4 = x4 - x5;
				const Zp z1 = y1 + y4, z4 = y1 - y4, z2 = y2 + y3, z3 = y2 - y3;
				x[k + 0 * n5] = z0 + u0;
				x[k + 1 * n5] = (z2 - u4) * w5[4 * k + 0];
				x[k + 2 * n5] = (z4 - u2) * w5[4 * k + 1];
				x[k + 3 * n5] = (z1 - u3) * w5[4 * k + 2];
				x[k + 4 * n5] = (z3 - u1) * w5[4 * k + 3];
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

	void nttinv(Zp * const x) const
	{
		const size_t n = _n, n5 = _n5;
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
			static const Zp K = Zp::primroot(5u), K1 = K + Zp(1u), K2 = K * K;
			static const Zp F2 = (K1 * K2).half(), F3 = K * K1, F4 = K * (K2 + Zp(1u)), F5 = (F3 + F4).half();

			// Radix-5
			for (size_t k = 0; k < n5; ++k)
			{
				const Zp u0 = x[k + 0 * n5], u4 = x[k + 1 * n5] * invw5[4 * k + 0], u3 = x[k + 2 * n5] * invw5[4 * k + 1];
				const Zp u2 = x[k + 3 * n5] * invw5[4 * k + 2], u1 = x[k + 4 * n5] * invw5[4 * k + 3];

				// 16 mul, 20 add
				// x[k + 0 * n5] = u0 + u1 + u2 + u3 + u4;
				// x[k + 1 * n5] = u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4;
				// x[k + 2 * n5] = u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4;
				// x[k + 3 * n5] = u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4;
				// x[k + 4 * n5] = u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4;

				// 4 mul, 20 add
				const Zp v1 = u1 + u4, v4 = u1 - u4, v2 = u2 + u3, v3 = u2 - u3;
				const Zp z0 = v1 + v2;
				const Zp x1 = u0, x2 = (v1 - v2) * F2, x3 = v3 * F3, x4 = v4 * F4, x5 = (v4 - v3) * F5;
				const Zp y1 = x1 + x2, y2 = x1 - x2, y3 = x3 + x5, y4 = x4 - x5;
				const Zp z1 = y1 + y4, z4 = y1 - y4, z2 = y2 + y3, z3 = y2 - y3;
				x[k + 0 * n5] = z0 + u0;
				x[k + 1 * n5] = z2 - u4;
				x[k + 2 * n5] = z4 - u2;
				x[k + 3 * n5] = z1 - u3;
				x[k + 4 * n5] = z3 - u1;
			}
		}
	}

	void mul1(Zp * const x, const Zp * const y) const
	{
		for (size_t k = 0, n = _n; k < n; ++k) x[k] *= y[k];
	}

	// Radix-2 NTT, mul, Radix-2 inverse NTT
	void mul2(Zp * const x, const Zp * const y) const
	{
		for (size_t k = 0, n = _n; k < n; k += 2)
		{
			const Zp u0 = x[k + 0], u1 = x[k + 1], up0 = y[k + 0], up1 = y[k + 1];
			const Zp v0 = u0 + u1, v1 = u0 - u1, vp0 = up0 + up1, vp1 = up0 - up1;
			const Zp s0 = v0 * vp0, s1 = v1 * vp1;
			x[k + 0] = s0 + s1; x[k + 1] = s0 - s1;
		}
	}

	void carry(Zp * const x) const
	{
		const size_t n = _n;
		const int * const digit_width = _digit_width;

		uint64_t c = 0;
		for (size_t k = 0; k < n; ++k)
		{
			x[k].set(digit_adc(x[k].get(), digit_width[k], c));
		} 

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
	Mersenne(const uint32_t q) : _n(transformsize(q)), _n5((_n % 5 == 0) ? _n / 5 : _n),
		_x(new Zp[3 * _n]),	// allocate 3 registers
		_w(new Zp[3 * _n5]), _invw(new Zp[3 * _n5]),
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
		const Zp nr2 = Zp(554u).pow((Zp::p() - 1) / 192 / n);
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

	virtual ~Mersenne()
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

	enum class Reg : uint32_t { R0 = 0, R1 = 1, R2 = 2 };

	void set(const Reg dst, const size_t a) const
	{
		const size_t n = _n;
		Zp * const x = &_x[uint32_t(dst) * n];

		x[0] = Zp(a);
		for (size_t k = 1; k < n; ++k) x[k] = Zp(0u);
	}

	void copy(const Reg dst, const Reg src) const
	{
		const size_t n = _n;
		const Zp * const x = &_x[uint32_t(src) * n];
		Zp * const y = &_x[uint32_t(dst) * n];

		for (size_t k = 0; k < n; ++k) y[k] = x[k];
	}

	bool is_equal(const uint32_t a) const
	{
		const size_t n = _n;
		const Zp * const x = _x;

		if (x[0].get() + (x[1].get() << _digit_width[0]) != a) return false;
		for (size_t k = 2; k < n; ++k) if (x[k].get() != 0) return false;
		return true;
	}

	bool is_equal(const Reg src1, const Reg src2) const
	{
		const size_t n = _n;
		const Zp * const x = &_x[uint32_t(src1) * n];
		const Zp * const y = &_x[uint32_t(src2) * n];

		for (size_t k = 0; k < n; ++k) if (y[k] != x[k]) return false;
		return true;
	}

	void square(const Reg src) const
	{
		Zp * const x = &_x[uint32_t(src) * _n];

		// weighted convolution
		mul1(x, _digit_weight);
		ntt(x);
		if (_even_exponent) mul1(x, x); else mul2(x, x);
		nttinv(x);
		mul1(x, _digit_invweight);

		// carry propagation
		carry(x);
	}

	// src is destroyed
	void mul(const Reg dst, const Reg src) const
	{
		Zp * const x = &_x[uint32_t(dst) * _n];
		Zp * const y = &_x[uint32_t(src) * _n];

		// weighted convolution
		mul1(y, _digit_weight);
		ntt(y);
		mul1(x, _digit_weight);
		ntt(x);
		if (_even_exponent) mul1(x, y); else mul2(x, y);
		nttinv(x);
		mul1(x, _digit_invweight);

		// carry propagation
		carry(x);
	}

	void error() const
	{
		_x[_n / 2] += Zp(1);
	}
};

int main()
{
	using Reg = Mersenne::Reg;

	size_t count = 0, count5 = 0;

	// 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, ...
	for (uint32_t p = 3; p <= 1207959503; p += 2)
	{
		bool isprime = true;
		for (uint32_t d = 3; p / d >= d; d += 2) if (p % d == 0) { isprime = false; break; }
		if (!isprime) continue;

		Mersenne m(p);

		++count;
		if (m.get_length() % 5 == 0) ++count5;

		// Gerbicz-Li error checking
		const uint32_t B_GL = std::max(uint32_t(std::sqrt(p)), 2u);

		// 3^{2^p}
		m.set(Reg::R0, 3u);	// result = 3
		m.set(Reg::R1, 1u);	// d(t) = 1
		for (uint32_t i = 0, j = p - 1; i < p; ++i, --j)
		{
			m.square(Reg::R0);
			if ((p == 1511) && (j == 0)) m.error();	// test Gerbicz-Li

			if ((j % B_GL == 0) && (j != 0))
			{
				m.copy(Reg::R2, Reg::R0);	// copy result
				m.mul(Reg::R1, Reg::R2);	// d(t + 1) = d(t) * result
			}
		}

		// 3-prp test: 3^{(2^p - 1) - 1} ?= 1  <=>  3^{2^p} ?= 9
		const bool is_prp = m.is_equal((p == 3) ? 9u % ((1u << p) - 1u) : 9u);

		// d(t + 1) = d(t) * result
		m.copy(Reg::R2, Reg::R1);
		m.mul(Reg::R2, Reg::R0);

		// The exponent of the residue is 2^(p mod B)
		// See: An Efficient Modular Exponentiation Proof Scheme, §2, Darren Li, Yves Gallot, https://arxiv.org/abs/2209.15623

		// d(t)^{2^B} * 3^{2^(p mod B)} = (3 * d(t)^{2^(B - p mod B)})^{2^(p mod B)}
		for (uint32_t i = 0; i < B_GL - p % B_GL; ++i) m.square(Reg::R1);
		m.set(Reg::R0, 3u);
		m.mul(Reg::R1, Reg::R0);
		for (uint32_t i = 0; i < p % B_GL; ++i) m.square(Reg::R1);

		if (!m.is_equal(Reg::R1, Reg::R2)) std::cout << p << ": Gerbicz-Li failed!" << std::endl;

		if (is_prp) std::cout << p << ": " << m.get_length() << " (radix-5: " << count5 << "/" << count << " = " << 100.0 * count5 / count <<  "%)" << std::endl;
	}

	return EXIT_SUCCESS;
}
