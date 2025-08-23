/*
Copyright 2025, Yves Gallot

mersenne2.cpp is free source code. You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <iostream>
#include <cstdint>

#define RIGHT_HANDED

// Z/{2^61 - 1}Z: the prime field of order p = 2^61 - 1
class Z61
{
private:
	static const uint64_t _p = (uint64_t(1) << 61) - 1;
	uint64_t _n;	// 0 <= n < p

	static uint64_t _add(const uint64_t a, const uint64_t b)
	{
		const uint64_t t = a + b;
		return t - ((t >= _p) ? _p : 0);
	}

	static uint64_t _sub(const uint64_t a, const uint64_t b)
	{
		const uint64_t t = a - b;
		return t + ((a < b) ? _p : 0);
	}

	static uint64_t _mul(const uint64_t a, const uint64_t b)
	{
		const __uint128_t t = a * __uint128_t(b);
		return _add(uint64_t(t) & _p, uint64_t(t >> 61));
	}

	static uint64_t _lshift(const uint64_t a, const uint8_t s)
	{
		const __uint128_t t = __uint128_t(a) << s;
		return _add(uint64_t(t) & _p, uint64_t(t >> 61));
	}

public:
	Z61() {}
	explicit Z61(const uint64_t n) : _n(n) {}

	uint64_t get() const { return _n; }

	bool operator!=(const Z61 & rhs) const { return (_n != rhs._n); }

	Z61 neg() const { return Z61((_n == 0) ? 0 : _p - _n); }

	Z61 operator+(const Z61 & rhs) const { return Z61(_add(_n, rhs._n)); }
	Z61 operator-(const Z61 & rhs) const { return Z61(_sub(_n, rhs._n)); }
	Z61 operator*(const Z61 & rhs) const { return Z61(_mul(_n, rhs._n)); }

	Z61 sqr() const { return Z61(_mul(_n, _n)); }

	Z61 lshift(const uint8_t s) const { const uint8_t s61 = s % 61; return (s61 != 0) ? Z61(_lshift(_n, s61)) : *this; }
	Z61 rshift(const uint8_t s) const { const uint8_t s61 = s % 61; return (s61 != 0) ? Z61(_lshift(_n, 61 - s61)) : *this; }
};

// GF((2^61 - 1)^2): the prime field of order p^2, p = 2^61 - 1
class GF61
{
private:
	Z61 _s0, _s1;
	// a primitive root of order 2^62 which is a root of (0, 1).
	static const uint64_t _h_order = uint64_t(1) << 62;
#ifdef RIGHT_HANDED
	static const uint64_t _h_0 = 481139922016222ull, _h_1 = 814659809902011ull;
#else
	static const uint64_t _h_0 = 264036120304204ull, _h_1 = 4677669021635377ull;
#endif

public:
	GF61() {}
	explicit GF61(const Z61 & s0, const Z61 & s1) : _s0(s0), _s1(s1) {}
	explicit GF61(const uint64_t n0, const uint64_t n1) : _s0(n0), _s1(n1) {}

	const Z61 & s0() const { return _s0; }
	const Z61 & s1() const { return _s1; }

	void set0(const uint64_t n0) { _s0 = Z61(n0); }
	void set1(const uint64_t n1) { _s1 = Z61(n1); }

	bool operator!=(const GF61 & rhs) const { return ((_s0 != rhs._s0) || (_s1 != rhs._s1)); }

	GF61 swap() const { return GF61(_s1, _s0); }
	GF61 conj() const { return GF61(_s0, _s1.neg()); }

	GF61 operator+(const GF61 & rhs) const { return GF61(_s0 + rhs._s0, _s1 + rhs._s1); }
	GF61 operator-(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, _s1 - rhs._s1); }
	GF61 addconj(const GF61 & rhs) const { return GF61(_s0 + rhs._s0, _s1 - rhs._s1); }
	GF61 subconj(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, _s1 + rhs._s1); }
	GF61 sub_conj(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, rhs._s1 - _s1); }
#ifdef RIGHT_HANDED
	GF61 addi(const GF61 & rhs) const { return GF61(_s0 + rhs._s1, _s1 - rhs._s0); }
	GF61 subi(const GF61 & rhs) const { return GF61(_s0 - rhs._s1, _s1 + rhs._s0); }
#else
	GF61 addi(const GF61 & rhs) const { return GF61(_s0 - rhs._s1, _s1 + rhs._s0); }
	GF61 subi(const GF61 & rhs) const { return GF61(_s0 + rhs._s1, _s1 - rhs._s0); }
#endif

GF61 sqr() const { const Z61 t = _s0 * _s1; return GF61(_s0.sqr() - _s1.sqr(), t + t); }
	GF61 mul(const GF61 & rhs) const { return GF61(_s0 * rhs._s0 - _s1 * rhs._s1, _s1 * rhs._s0 + _s0 * rhs._s1); }
	GF61 mulconj(const GF61 & rhs) const { return GF61(_s0 * rhs._s0 + _s1 * rhs._s1, _s1 * rhs._s0 - _s0 * rhs._s1); }

	GF61 lshift(const uint8_t ls0, const uint8_t ls1) const { return GF61(_s0.lshift(ls0), _s1.lshift(ls1)); }
	GF61 rshift(const uint8_t rs0, const uint8_t rs1) const { return GF61(_s0.rshift(rs0), _s1.rshift(rs1)); }

	GF61 pow(const uint64_t e) const
	{
		if (e == 0) return GF61(1u, 0u);
		GF61 r = GF61(1u, 0u), y = *this;
		for (uint64_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r = r.mul(y); y = y.sqr(); }
		return r.mul(y);
	}

	static const GF61 root_nth(const size_t n) { return GF61(Z61(_h_0), Z61(_h_1)).pow(_h_order / n); }
	static uint8_t log2_root_two(const size_t n) { return uint8_t(((uint64_t(1) << 60) / n) % 61); }
};

class mersenne
{
private:
	const uint8_t _ln;
	const size_t _n;
	const bool _even;
	GF61 * const _z;
	GF61 * const _w;
	uint8_t * const _w_ib;
	uint8_t * const _digit_width;

private:
	static constexpr bool is_even(const size_t n)
	{
		size_t m = n; for (; m > 1; m /= 4);
		return (m == 1);
	}

	// Make sure the transform is long enough so that each digit cannot overflow 2^61 - 1 after the convolution.
	static constexpr uint8_t transformsize(const uint32_t q)
	{
		uint8_t ln = 2; uint32_t w = 0;
		do
		{
			++ln;
			// digit-width is w or w + 1
			w = q >> ln;
		// the condition is n * (2^{w + 1} - 1)^2 < 2^61 - 1
		} while (ln + 2 * (w + 1) >= 61);

		return ln;
	}

	void swap() const
	{
		GF61 * const z = _z;
		for (size_t i = 0, n_2 = _n / 2; i < n_2; ++i) z[i] = z[i].swap();
	}

	void conjugate() const
	{
		GF61 * const z = _z;
		for (size_t i = 0, n_2 = _n / 2; i < n_2; ++i) z[i] = z[i].conj();
	}

	// Bit-reversed permutation
	void scramble() const
	{
		GF61 * const z = _z;

		for (size_t i = 0, j = 0, n_2 = _n / 2; i < n_2 - 1; ++i)
		{
			if (i < j) std::swap(z[i], z[j]);
			size_t k = n_2 / 2; while (k <= j) { j -= k; k /= 2; }
			j += k;
		}
	}

	// Radix-2, Decimation-In-Frequency
	void DIF2(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 wj = w[j * s];

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 2 * m * i + j;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m];
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = (u0 - u1).mulconj(wj);
			}
		}
	}

	// Radix-2, Decimation-In-Time
	void DIT2(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 wj = w[j * s];

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 2 * m * i + j;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mulconj(wj);
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = u0 - u1;
			}
		}
	}

	// Inverse radix-2, Decimation-In-Frequency
	void InverseDIF2(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 wj = w[j * s];

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 2 * m * i + j;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m];
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = (u0 - u1).mul(wj);
			}
		}
	}

	// Inverse Radix-2, Decimation-In-Time
	void InverseDIT2(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 wj = w[j * s];

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 2 * m * i + j;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mul(wj);
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = u0 - u1;
			}
		}
	}

	// Radix-4, Decimation-In-Frequency
	void DIF4(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 w1 = w[j * s], w2 = w[2 * j * s], w3 = w1.mul(w2);

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 4 * m * i + j;

				const GF61 u0 = z[k + 0 * m], u2 = z[k + 2 * m], u1 = z[k + 1 * m], u3 = z[k + 3 * m];
				const GF61 v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = u1 - u3;
				z[k + 0 * m] = v0 + v1;
				z[k + 1 * m] = (v0 - v1).mulconj(w2);
				z[k + 2 * m] = v2.subi(v3).mulconj(w1);
				z[k + 3 * m] = v2.addi(v3).mulconj(w3);
			}
		}
	}

	// Radix-4, Decimation-In-Time
	void DIT4(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 w1 = w[j * s], w2 = w[2 * j * s], w3 = w1.mul(w2);

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 4 * m * i + j;

				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mulconj(w2), u2 = z[k + 2 * m].mulconj(w1), u3 = z[k + 3 * m].mulconj(w3);
				const GF61 v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = u3 - u2;
				z[k + 0 * m] = v0 + v2;
				z[k + 2 * m] = v0 - v2;
				z[k + 1 * m] = v1.addi(v3);
				z[k + 3 * m] = v1.subi(v3);
			}
		}
	}

	// Inverse radix-4, Decimation-In-Frequency
	void InverseDIF4(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 w1 = w[j * s], w2 = w[2 * j * s], w3 = w1.mul(w2);

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 4 * m * i + j;

				const GF61 u0 = z[k + 0 * m], u2 = z[k + 2 * m], u1 = z[k + 1 * m], u3 = z[k + 3 * m];
				const GF61 v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = u1 - u3;
				z[k + 0 * m] = v0 + v1;
				z[k + 1 * m] = (v0 - v1).mul(w2);
				z[k + 2 * m] = v2.addi(v3).mul(w1);
				z[k + 3 * m] = v2.subi(v3).mul(w3);
			}
		}
	}

	// Inverse radix-4, Decimation-In-Time
	void InverseDIT4(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < m; ++j)
		{
			const GF61 w1 = w[j * s], w2 = w[2 * j * s], w3 = w1.mul(w2);

			for (size_t i = 0; i < s; ++i)
			{
				const size_t k = 4 * m * i + j;

				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mul(w2), u2 = z[k + 2 * m].mul(w1), u3 = z[k + 3 * m].mul(w3);
				const GF61 v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = u3 - u2;
				z[k + 0 * m] = v0 + v2;
				z[k + 2 * m] = v0 - v2;
				z[k + 1 * m] = v1.subi(v3);
				z[k + 3 * m] = v1.addi(v3);
			}
		}
	}

	// Decimation-In-Frequency
	void NTT_DIF() const
	{
		for (size_t m = _n / 8, s = 1; m >= 1; m /= 4, s *= 4) DIF4(m, s);
		if (_even) DIF2(1, _n / 4);
		scramble();
	}

	// Decimation-In-Time
	void NTT_DIT() const
	{
		scramble();
		for (size_t m = 1, s = _n / 8; s >= 1; m *= 4, s /= 4) DIT4(m, s);
		if (_even) DIT2(_n / 4, 1);
	}

	// Inverse Decimation-In-Frequency
	void InverseNTT_DIF() const
	{
		for (size_t m = _n / 8, s = 1; m >= 1; m /= 4, s *= 4) InverseDIF4(m, s);
		if (_even) InverseDIF2(1, _n / 4);
		scramble();
	}

	// Inverse Decimation-In-Time
	void InverseNTT_DIT() const
	{
		scramble();
		for (size_t m = 1, s = _n / 8; s >= 1; m *= 4, s /= 4) InverseDIT4(m, s);
		if (_even) InverseDIT2(_n / 4, 1);
	}

	// IBDWT: weighted digits
	void weight() const
	{
		GF61 * const z = _z;
		const uint8_t * const w_ib = _w_ib;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
			z[k] = z[k].lshift(w_ib[2 * k + 0], w_ib[2 * k + 1]);
		}
	}

	// IBDWT: restore the unweighted digits. sqr and backward transform must be divided by 2n
	void unweight_norm() const
	{
		const uint8_t ln1 = _ln + 1;
		GF61 * const z = _z;
		const uint8_t * const w_ib = _w_ib;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
			z[k] = z[k].rshift(w_ib[2 * k + 0] + ln1, w_ib[2 * k + 1] + ln1);
		}
	}

	// Input of the transform is 'real' (if alpha is the symbolic square root of p - 1 and the elements of GF(p^2) are a + b*alpha then b = 0).
	// A length-n/2 transform is computed onto z(k) = x(2*k + 0) + x(2*k + 1)*alpha.
	// See Henrik V. Sorensen, Douglas L. Jones, Michael T. Heideman, C. Sidney Burrus, "Real-Valued Fast Fourier Transform Algorithms",
	// in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 35, no. 6, pp. 849-863, June 1987.
	// Output is recombined to produce the half-length transform (full length is not needed because of Hermitian symmetry).
	void sqr() const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;
		const size_t n_4 = _n / 4;

		const GF61 u = z[0] + z[0];
		z[0] = GF61(u.s0().sqr() + u.s1().sqr(), u.s0() * (u.s1() + u.s1()));
		z[n_4] = (z[n_4] + z[n_4]).sqr();

		for (size_t j = 1, mj = _n / 2 - 1; j < n_4; ++j, --mj)
		{
			const GF61 w2 = w[j];
			const GF61 zj = z[j], zmj = z[mj];
			const GF61 u0 = zj.addconj(zmj), u1 = zj.subconj(zmj);
#ifdef RIGHT_HANDED
			const GF61 v0 = u0.sqr() - u1.sqr().mulconj(w2), v1 = u0.mul(u1 + u1);
#else
			const GF61 v0 = u0.sqr() - u1.sqr().mul(w2), v1 = u0.mul(u1 + u1);
#endif
			z[j] = v0 + v1; z[mj] = v0.sub_conj(v1);
		}
	}

	// Add a carry to the number and return the carry of the first digit_width bits
	static constexpr uint64_t digit_adc(const uint64_t lhs, const uint8_t digit_width, uint64_t & carry)
	{
		const uint64_t s = lhs + carry;
		carry = s >> digit_width;
		return s & ((uint64_t(1) << digit_width) - 1);
	}

	// Adjust the digits to the digit representation
	void carry() const
	{
		const size_t n_2 = _n / 2;
		GF61 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		uint64_t c = 0;
		for (size_t k = 0; k < n_2; ++k)
		{
			const uint64_t n0 = digit_adc(z[k].s0().get(), digit_width[2 * k + 0], c);
			const uint64_t n1 = digit_adc(z[k].s1().get(), digit_width[2 * k + 1], c);
			z[k] = GF61(n0, n1);
		}

		while (c != 0)
		{
			for (size_t k = 0; k < n_2; ++k)
			{
				const uint64_t n0 = digit_adc(z[k].s0().get(), digit_width[2 * k + 0], c);
				const uint64_t n1 = digit_adc(z[k].s1().get(), digit_width[2 * k + 1], c);
				z[k] = GF61(n0, n1);
				if (c == 0) break;
			}
		}
	}

public:
	mersenne(const uint32_t q) : _ln(transformsize(q)), _n(size_t(1) << _ln), _even(is_even(_n)),
		_z(new GF61[_n / 2]), _w(new GF61[_n / 4]), _w_ib(new uint8_t[_n]), _digit_width(new uint8_t[_n])
	{
		const uint8_t ln = _ln;
		const size_t n = _n;

		GF61 * const w = _w;
		const GF61 r = GF61::root_nth(n / 2);

		// Radix-2 twiddle factors
		GF61 rj = GF61(1, 0);
		for (size_t j = 0; j < n / 4; ++j) { w[j] = rj; rj = rj.mul(r); }

		// IBDWT weights: x^q - 1 => x^n - 1
		// See Richard Crandall, Barry Fagin, "Discrete weighted transforms and large-integer arithmetic", Math. Comp. 62 (1994), 305-324.

		// Weights are power of two. Store log_2(weight).
		uint8_t * const w_ib = _w_ib;
		uint8_t * const digit_width = _digit_width;

		const uint8_t lr2_61 = GF61::log2_root_two(n);		// n-th root of two

		w_ib[0] = 0;

		const uint8_t q_n = uint8_t(q >> ln);
		// std::cout << q << ": length = " << n << ", digit-width = " << int(q_n) << "/" << int(q_n + 1) << std::endl;

		uint32_t o = 0;
		for (size_t j = 0; j <= n; ++j)
		{
			const uint64_t qj = q * uint64_t(j);
			// ceil(a / b) = floor((a - 1) / b) + 1
			const uint32_t ceil_qj_n = uint32_t(((qj - 1) >> ln) + 1);

			// bit position for digit[i] is ceil(qj / n)
			if (j > 0)
			{
				const uint8_t c = uint8_t(ceil_qj_n - o);
				if ((c != q_n) && (c != q_n + 1)) throw;
				digit_width[j - 1] = c;

				// weight is 2^[ceil(qj / n) - qj/n]
				if (j < n)
				{
					// e = (ceil(qj / n).n - qj) / n
					// qj = k.n => e = 0
					// qj = k.n + r, r > 0 => ((k + 1).n - (k.n + r)) / n = (n - r) / n
					const uint32_t r = uint32_t(qj & (n - 1));
					w_ib[j] = uint8_t((lr2_61 * (n - r)) % 61);
				}
			}

			o = ceil_qj_n;
		}
	}

	virtual ~mersenne()
	{
		delete[] _z;
		delete[] _w;
		delete[] _w_ib;
		delete[] _digit_width;
	}

	size_t get_n() const { return _n; }

	void init(const uint64_t a) const
	{
		GF61 * const z = _z;

		z[0] = GF61(a, 0);
		for (size_t k = 1, n_2 = _n / 2; k < n_2; ++k) z[k] = GF61(0, 0);
	}

	void square() const
	{
		// weighted convolution, radix-4 transform
		weight();
#ifdef RIGHT_HANDED
		NTT_DIF();
		// NTT_DIT();
#else
		InverseNTT_DIF();
		// InverseNTT_DIT();
#endif
		sqr();
#ifdef RIGHT_HANDED
		// InverseNTT_DIF();
		// InverseNTT_DIT();
		conjugate(); NTT_DIF(); conjugate();
		// swap(); NTT_DIF(); swap();
#else
		// NTT_DIF();
		// NTT_DIT();
		conjugate(); InverseNTT_DIF(); conjugate();
		// swap(); InverseNTT_DIF(); swap();
#endif
		unweight_norm();

		// carry propagation
		carry();
	}

	// Subtract a carry to the number and return the carry if borrowing
	static constexpr uint64_t digit_sbc(const uint64_t lhs, const uint8_t digit_width, uint32_t & carry)
	{
		const bool borrow = (lhs < carry);
		const uint64_t r = lhs - carry + (borrow ? (uint64_t(1) << digit_width) : 0);
		carry = borrow ? 1 : 0;
		return r;
	}

	void sub(const uint32_t a) const
	{
		GF61 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		uint32_t c = a;
		while (c != 0)
		{
			for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
			{
				const uint64_t n0 = digit_sbc(z[k].s0().get(), digit_width[2 * k + 0], c);
				const uint64_t n1 = digit_sbc(z[k].s1().get(), digit_width[2 * k + 1], c);
				z[k] = GF61(n0, n1);
				if (c == 0) break;
			}
		}
	}

	bool is_zero() const
	{
		const GF61 * const z = _z;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
			if (z[k].s0().get() != 0u) return false;
			if (z[k].s1().get() != 0u) return false;
		}
		return true;
	}

	bool is_Mp() const
	{
		const GF61 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
			const uint32_t n0 = z[k].s0().get();
			if (n0 != (uint64_t(1) << digit_width[2 * k + 0]) - 1) return false;
			const uint32_t n1 = z[k].s1().get();
			if (n1 != (uint64_t(1) << digit_width[2 * k + 1]) - 1) return false;
		}

		return true;
	}
};

int main()
{
	// 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, ...
	for (uint32_t p = 3; p <= 4294967291; p += 2)
	{
		bool isprime = true;
		for (uint32_t d = 3; p / d >= d; d += 2) if (p % d == 0) { isprime = false; break; }
		if (!isprime) continue;

		mersenne m(p);
		// std::cout << p << ", " << m.get_n() << std::endl;

		// Lucas–Lehmer primality test
		m.init(4);
		for (uint32_t i = 0; i < p - 2; ++i)
		{
			m.square();
			m.sub(2);
		}

		// IBDWT is modulo 2^p then 0 (mod p) is 0 or 2^p - 1.
		if (m.is_zero() || m.is_Mp()) std::cout << p << std::endl;
	}

	return EXIT_SUCCESS;
}
