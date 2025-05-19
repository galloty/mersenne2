/*
Copyright 2025, Yves Gallot

mersenne2.cpp is free source code. You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <iostream>
#include <cstdint>

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
		const uint64_t lo = uint64_t(t), hi = uint64_t(t >> 64);
		const uint64_t lo61 = lo & _p, hi61 = (lo >> 61) | (hi << 3);
		return _add(lo61, hi61);
	}

	static uint64_t _lshift(const uint64_t a, const uint8_t s)
	{
		const uint64_t lo = a << s, hi = a >> (64 - s);
		const uint64_t lo61 = lo & _p, hi61 = (lo >> 61) | (hi << 3);
		return _add(lo61, hi61);
	}

public:
	Z61() {}
	explicit Z61(const uint64_t n) : _n(n) {}

	uint64_t get() const { return _n; }

	bool operator!=(const Z61 & rhs) const { return (_n != rhs._n); }

	// Z61 neg() const { return Z61((_n == 0) ? 0 : _p - _n); }
	// Z61 half() const { return Z61(((_n % 2 == 0) ? _n : (_n + _p)) / 2); }

	Z61 operator+(const Z61 & rhs) const { return Z61(_add(_n, rhs._n)); }
	Z61 operator-(const Z61 & rhs) const { return Z61(_sub(_n, rhs._n)); }
	Z61 operator*(const Z61 & rhs) const { return Z61(_mul(_n, rhs._n)); }

	Z61 sqr() const { return Z61(_mul(_n, _n)); }

	Z61 lshift(const uint8_t s) const { const uint8_t s61 = s % 61; return (s61 != 0) ? Z61(_lshift(_n, s61)) : *this; }
	Z61 rshift(const uint8_t s) const { const uint8_t s61 = s % 61; return (s61 != 0) ? Z61(_lshift(_n, 61 - s61)) : *this; }

	Z61 pow(const uint64_t e) const
	{
		if (e == 0) return Z61(1u);
		Z61 r = Z61(1u), y = *this;
		for (uint64_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r = r * y; y = y.sqr(); }
		return r * y;
	}

	Z61 invert() const { return Z61(pow(_p - 2)); }
};

// GF((2^61 - 1)^2): the prime field of order p^2, p = 2^61 - 1
class GF61
{
private:
	Z61 _s0, _s1;
	// a primitive root of order 2^62 which is a root of (0, 1).
	static const uint64_t _h_order = uint64_t(1) << 62;
	static const uint64_t _h_0 = 264036120304204ull, _h_1 = 4677669021635377ull;

public:
	GF61() {}
	explicit GF61(const Z61 & s0, const Z61 & s1) : _s0(s0), _s1(s1) {}
	explicit GF61(const uint64_t n0, const uint64_t n1) : _s0(n0), _s1(n1) {}

	const Z61 & s0() const { return _s0; }
	const Z61 & s1() const { return _s1; }

	void set0(const uint64_t l0) { _s0 = Z61(l0); }
	void set1(const uint64_t l1) { _s1 = Z61(l1); }

	bool operator!=(const GF61 & rhs) const { return ((_s0 != rhs._s0) || (_s1 != rhs._s1)); }

	// GF61 conj() const { return GF61(_s0, _s1.neg()); }
	// GF61 muli() const { return GF61(-_s1, _s0); }
	// GF61 half() const { return GF61(_s0.half(), _s1.half()); }

	GF61 operator+(const GF61 & rhs) const { return GF61(_s0 + rhs._s0, _s1 + rhs._s1); }
	GF61 operator-(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, _s1 - rhs._s1); }
	GF61 addconj(const GF61 & rhs) const { return GF61(_s0 + rhs._s0, _s1 - rhs._s1); }
	GF61 subconj(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, _s1 + rhs._s1); }
	GF61 addi(const GF61 & rhs) const { return GF61(_s0 - rhs._s1, _s1 + rhs._s0); }
	GF61 subi(const GF61 & rhs) const { return GF61(_s0 + rhs._s1, _s1 - rhs._s0); }
	GF61 subi_conj(const GF61 & rhs) const { return GF61(_s0 + rhs._s1, rhs._s0 - _s1); }

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

	static const GF61 root_one(const size_t n) { return GF61(Z61(_h_0), Z61(_h_1)).pow(_h_order / n); }
	static uint8_t log2_root_two(const size_t n) { return uint8_t(((uint64_t(1) << 60) / n) % 61); }

	// Add a carry to the number and return the carry of the first digit_width bits
	static constexpr uint64_t digit_adc(const uint64_t lhs, const uint8_t digit_width, uint64_t & carry)
	{
		const uint64_t s = lhs + carry;
		const uint64_t c = (s < lhs) ? 1 : 0;
		carry = (s >> digit_width) + (c << (64 - digit_width));
		return s & ((uint64_t(1) << digit_width) - 1);
	}

	// Subtract a carry to the number and return the carry of borrowing
	static constexpr uint64_t digit_sbc(const uint64_t lhs, const uint8_t digit_width, uint64_t & carry)
	{
		const bool borrow = (lhs < carry);
		const uint64_t r = lhs - carry + (borrow ? (uint64_t(1) << digit_width) : 0);
		carry = borrow ? 1 : 0;
		return r;
	}

	// Add a carry
	GF61 adc(const uint8_t digit_width0, const uint8_t digit_width1, uint64_t & carry) const
	{
		const uint64_t d0 = digit_adc(_s0.get(), digit_width0, carry);
		const uint64_t d1 = digit_adc(_s1.get(), digit_width1, carry);
		return GF61(d0, d1);
	}

	// Subtract a carry
	GF61 sbc(const uint8_t digit_width0, const uint8_t digit_width1, uint64_t & carry) const
	{
		const uint64_t d0 = digit_sbc(_s0.get(), digit_width0, carry);
		if (carry == 0) return GF61(d0, _s1.get());
		const uint64_t d1 = digit_sbc(_s1.get(), digit_width1, carry);
		return GF61(d0, d1);
	}
};

class mersenne
{
private:
	const uint8_t _ln;
	const size_t _n;
	GF61 * const _z;
	GF61 * const _w;
	uint8_t * const _w_ib;
	uint8_t * const _digit_width;

private:
	// bit-reversal permutation of index i for a sequence of n items
	static constexpr size_t bitrev(const size_t i, const size_t n)
	{
		size_t r = 0;
		for (size_t k = n, j = i; k > 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
		return r;
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
		// the condition is n * (2^{w + 1} - 1)^2 < 2^61 - 1)
		} while (ln + 2 * (w + 1) >= 61);

		return ln;
	}

	// Bruun's recursive polynomial factorization
	// See G. Bruun, "z-transform DFT filters and FFT's," in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 26, no. 1, pp. 56-63, February 1978.

	/*void forward2(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < s; ++j)
		{
			const GF61 wj = w[s + j];
			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 2 * m * j + i;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mul(wj);
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = u0 - u1;
			}
		}
	}*/

	void forward4(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;
		const GF61 * const wr4 = &_w[_n / 2];

		for (size_t j = 0; j < s; ++j)
		{
			const GF61 w1 = w[s + j], w2 = w[2 * (s + j)], w3 = wr4[s + j];

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 4 * m * j + i;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mul(w2), u2 = z[k + 2 * m].mul(w1), u3 = z[k + 3 * m].mul(w3);
				const GF61 v0 = u0 + u2, v1 = u1 + u3, v2 = u0 - u2, v3 = u1 - u3;
				z[k + 0 * m] = v0 + v1; z[k + 1 * m] = v0 - v1;
				z[k + 2 * m] = v2.addi(v3); z[k + 3 * m] = v2.subi(v3);
			}
		}
	}

	// Inverse transform of Bruun's method

	/*void backward2(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < s; ++j)
		{
			const GF61 wj = w[s + j];
			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 2 * m * j + i;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m];
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = (u0 - u1).mulconj(wj);
			}
		}
	}*/

	void backward4(const size_t m, const size_t s) const
	{
		GF61 * const z = _z;
		const GF61 * const w = _w;
		const GF61 * const wr4 = &_w[_n / 2];

		for (size_t j = 0; j < s; ++j)
		{
			const GF61 w1 = w[s + j], w2 = w[2 * (s + j)], w3 = wr4[s + j];

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 4 * m * j + i;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m], u2 = z[k + 2 * m], u3 = z[k + 3 * m];
				const GF61 v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = u3 - u2;
				z[k + 0 * m] = v0 + v2; z[k + 2 * m] = (v0 - v2).mulconj(w1);
				z[k + 1 * m] = v1.addi(v3).mulconj(w2); z[k + 3 * m] = v1.subi(v3).mulconj(w3);
			}
		}
	}

	// IBDWT: weighted digits
	void weight() const
	{
		GF61 * const z = _z;
		const uint8_t * const w_ib = _w_ib;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k) z[k] = z[k].lshift(w_ib[2 * k + 0], w_ib[2 * k + 1]);
	}

	// IBDWT: restore the unweighted digits. sqr and backward transform must be divided by 4n
	void unweight_norm() const
	{
		const uint8_t ln = _ln;
		GF61 * const z = _z;
		const uint8_t * const w_ib = _w_ib;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k) z[k] = z[k].rshift(w_ib[2 * k + 0] + ln + 2, w_ib[2 * k + 1] + ln + 2);
	}

	// Input of the transform is 'real' (if alpha is the symbolic square root of p - 1 and the elements of GF(p^2) are a + b*alpha then b = 0).
	// A length-n/2 transform is computed onto z(k) = x(2*k + 0) + x(2*k + 1)*alpha.
	// See Henrik V. Sorensen, Douglas L. Jones, Michael T. Heideman, C. Sidney Burrus, "Real-Valued Fast Fourier Transform Algorithms",
	// in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 35, no. 6, pp. 849-863, June 1987.
	// Output is recombined to produce the half-length transform (full length is not needed because of Hermitian symmetry).
	void sqr() const
	{
		const size_t n_2 = _n / 2;
		GF61 * const z = _z;
		const GF61 * const w = &_w[_n];

		const GF61 Z0 = z[0] + z[0];
		const Z61 X0 = Z0.s0() + Z0.s1(), Xm0 = Z0.s0() - Z0.s1();
		const Z61 X0sq = X0.sqr(), Xm0sq = Xm0.sqr();
		z[0] = GF61(X0sq + Xm0sq, X0sq - Xm0sq);

		const GF61 zn_4 = z[1] + z[1];
		const GF61 Zn_4 = zn_4.sqr();
		z[1] = Zn_4 + Zn_4;

		for (size_t k = 2; k < n_2; k += 2)
		{
			const GF61 wk = w[k / 2];

			// const size_t j = bitrev(k, n_2), mk = bitrev(n_2 - j, n_2);
			const size_t mk = (size_t(3) << (63 - __builtin_clzll((unsigned long long)k))) - k - 1;
			const GF61 zk = z[k], zmk = z[mk];

			const GF61 Zek = zk.addconj(zmk), Zok = zk.subconj(zmk).mul(wk);
			const GF61 Zk = Zek.subi(Zok), Zmk = Zek.addi(Zok);

			const GF61 Zk2 = Zk.sqr(), Zmk2 = Zmk.sqr();

			const GF61 Zek2 = Zk2 + Zmk2, Zok2 = (Zk2 - Zmk2).mulconj(wk);
			const GF61 zk2 =  Zek2.addi(Zok2), zmk2 = Zek2.subi_conj(Zok2);

			z[k] = zk2; z[mk] = zmk2;
		}
	}

	// forward2, sqr, backward2
	void sqr2() const
	{
		const size_t n_4 = _n / 4;
		GF61 * const z = _z;
		const GF61 * const w = _w;

		for (size_t j = 0; j < n_4; ++j)
		{
			const GF61 u0 = z[2 * j + 0], u1 = z[2 * j + 1].mul(w[n_4 + j]);
			z[2 * j + 0] = u0 + u1;
			z[2 * j + 1] = u0 - u1;
		}

		sqr();

		for (size_t j = 0; j < n_4; ++j)
		{
			const GF61 u0 = z[2 * j + 0], u1 = z[2 * j + 1];
			z[2 * j + 0] = u0 + u1;
			z[2 * j + 1] = (u0 - u1).mulconj(w[n_4 + j]);
		}
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
			z[k] = z[k].adc(digit_width[2 * k + 0], digit_width[2 * k + 1], c);
		}

		while (c != 0)
		{
			for (size_t k = 0; k < n_2; ++k)
			{
				z[k] = z[k].adc(digit_width[2 * k + 0], digit_width[2 * k + 1], c);
				if (c == 0) break;
			}
		}
	}

public:
	mersenne(const uint32_t q) : _ln(transformsize(q)), _n(size_t(1) << _ln),
		_z(new GF61[_n]), _w(new GF61[5 * _n / 4]), _w_ib(new uint8_t[_n]), _digit_width(new uint8_t[_n])
	{
		const size_t n = _n;

		// radix-2 twiddle factors
		GF61 * const w = _w;
		for (size_t s = 1; s <= n / 4; s *= 2)
		{
			const GF61 r_s = GF61::root_one(2 * s);
			for (size_t j = 0; j < s; ++j) w[s + j] = r_s.pow(bitrev(j, s));
		}

		// radix-4 twiddle factors
		for (size_t s = 1; s <= n / 4; s *= 2)
		{
			for (size_t j = 0; j < s; ++j) w[n / 2 + s + j] = w[s + j].mul(w[2 * (s + j)]);
		}

		// n values in GF(p) and a transform of length n/2 in GF(p^2) 
		const GF61 r_n = GF61::root_one(n);
		for (size_t j = 0, s = n / 4; j < s; ++j) w[n + j] = r_n.pow(bitrev(j, s));

		// IBDWT weights: x^q - 1 => x^n - 1
		// See Richard Crandall, Barry Fagin, "Discrete weighted transforms and large-integer arithmetic", Math. Comp. 62 (1994), 305-324.

		// Weights are power of two. Store log_2(weight).
		uint8_t * const w_ib = _w_ib;
		uint8_t * const digit_width = _digit_width;

		const uint8_t lr2 = GF61::log2_root_two(n);		// n-th root of two

		w_ib[0] = 0;

		const uint8_t q_n = uint8_t(q / n);
		uint32_t o = 0;
		for (size_t j = 0; j <= n; ++j)
		{
			const uint64_t qj = q * uint64_t(j);
			// ceil(a / b) = floor((a - 1) / b) + 1
			const uint32_t ceil_qj_n = uint32_t((qj - 1) / n + 1);

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
					const uint32_t r = uint32_t(qj % n);
					w_ib[j] = uint8_t((lr2 * (n - r)) % 61);
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

	void init(const uint64_t a) const
	{
		GF61 * const z = _z;

		z[0] = GF61(a, 0u);
		for (size_t k = 1, n = _n; k < n / 2; ++k) z[k] = GF61(0u, 0u);
	}

	void square() const
	{
		// weighted convolution, radix-4 transform
		weight();
		size_t m = _n / 4, s = 1;
		for (; m > 1; m /= 4, s *= 4) forward4(m/2, s);
		if (m == 1) sqr2(); else sqr();
		for (m = (m == 1) ? 4 : 2, s /= 4; s >= 1; m *= 4, s /= 4) backward4(m/2, s);
		unweight_norm();

		// carry propagation
		carry();
	}

	void sub(const uint64_t a) const
	{
		GF61 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		uint64_t c = a;
		while (c != 0)
		{
			for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
			{
				z[k] = z[k].sbc(digit_width[2 * k + 0], digit_width[2 * k + 1], c);
				if (c == 0) break;
			}
		}
	}

	bool is_zero() const
	{
		GF61 * const z = _z;
		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k) if (z[k] != GF61(0u, 0u)) return false;
		return true;
	}

	bool is_Mp() const
	{
		GF61 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
 			if (z[k] != GF61((uint64_t(1) << digit_width[2 * k + 0]) - 1, (uint64_t(1) << digit_width[2 * k + 1]) - 1)) return false;
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
