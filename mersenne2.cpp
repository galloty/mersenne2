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

	// Z61 neg() const { return Z61((_n == 0) ? 0 : _p - _n); }
	// Z61 half() const { return Z61(((_n % 2 == 0) ? _n : (_n + _p)) / 2); }

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
	static const uint64_t _h_0 = 264036120304204ull, _h_1 = 4677669021635377ull;

public:
	GF61() {}
	explicit GF61(const Z61 & s0, const Z61 & s1) : _s0(s0), _s1(s1) {}
	explicit GF61(const uint64_t n0, const uint64_t n1) : _s0(n0), _s1(n1) {}

	const Z61 & s0() const { return _s0; }
	const Z61 & s1() const { return _s1; }

	void set0(const uint64_t n0) { _s0 = Z61(n0); }
	void set1(const uint64_t n1) { _s1 = Z61(n1); }

	bool operator!=(const GF61 & rhs) const { return ((_s0 != rhs._s0) || (_s1 != rhs._s1)); }

	// GF61 conj() const { return GF61(_s0, _s1.neg()); }
	// GF61 muli() const { return GF61(_s1.neg(), _s0); }
	// GF61 half() const { return GF61(_s0.half(), _s1.half()); }

	GF61 operator+(const GF61 & rhs) const { return GF61(_s0 + rhs._s0, _s1 + rhs._s1); }
	GF61 operator-(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, _s1 - rhs._s1); }
	GF61 addconj(const GF61 & rhs) const { return GF61(_s0 + rhs._s0, _s1 - rhs._s1); }
	GF61 subconj(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, _s1 + rhs._s1); }
	GF61 sub_conj(const GF61 & rhs) const { return GF61(_s0 - rhs._s0, rhs._s1 - _s1); }
	GF61 addi(const GF61 & rhs) const { return GF61(_s0 - rhs._s1, _s1 + rhs._s0); }
	GF61 subi(const GF61 & rhs) const { return GF61(_s0 + rhs._s1, _s1 - rhs._s0); }

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

// Z/{2^31 - 1}Z: the prime field of order p = 2^31 - 1
class Z31
{
private:
	static const uint32_t _p = (uint32_t(1) << 31) - 1;
	uint32_t _n;	// 0 <= n < p

	static uint32_t _add(const uint32_t a, const uint32_t b)
	{
		const uint32_t t = a + b;
		return t - ((t >= _p) ? _p : 0);
	}

	static uint32_t _sub(const uint32_t a, const uint32_t b)
	{
		const uint32_t t = a - b;
		return t + ((a < b) ? _p : 0);
	}

	static uint32_t _mul(const uint32_t a, const uint32_t b)
	{
		const uint64_t t = a * uint64_t(b);
		return _add(uint32_t(t) & _p, uint32_t(t >> 31));
	}

	static uint32_t _lshift(const uint32_t a, const uint8_t s)
	{
		const uint64_t t = uint64_t(a) << s;
		return _add(uint32_t(t) & _p, uint32_t(t >> 31));
	}

public:
	Z31() {}
	explicit Z31(const uint32_t n) : _n(n) {}
	explicit Z31(const uint64_t n) : _n(n % _p) {}

	uint32_t get() const { return _n; }

	bool operator!=(const Z31 & rhs) const { return (_n != rhs._n); }

	// Z31 neg() const { return Z31((_n == 0) ? 0 : _p - _n); }
	// Z31 half() const { return Z31(((_n % 2 == 0) ? _n : (_n + _p)) / 2); }

	Z31 operator+(const Z31 & rhs) const { return Z31(_add(_n, rhs._n)); }
	Z31 operator-(const Z31 & rhs) const { return Z31(_sub(_n, rhs._n)); }
	Z31 operator*(const Z31 & rhs) const { return Z31(_mul(_n, rhs._n)); }

	Z31 sqr() const { return Z31(_mul(_n, _n)); }

	Z31 lshift(const uint8_t s) const { const uint8_t s31 = s % 31; return (s31 != 0) ? Z31(_lshift(_n, s31)) : *this; }
	Z31 rshift(const uint8_t s) const { const uint8_t s31 = s % 31; return (s31 != 0) ? Z31(_lshift(_n, 31 - s31)) : *this; }
};

// GF((2^31 - 1)^2): the prime field of order p^2, p = 2^31 - 1
class GF31
{
private:
	Z31 _s0, _s1;
	// a primitive root of order 2^32 which is a root of (0, 1).
	static const uint64_t _h_order = uint64_t(1) << 32;
	static const uint32_t _h_0 = 7735u, _h_1 = 748621u;

public:
	GF31() {}
	explicit GF31(const Z31 & s0, const Z31 & s1) : _s0(s0), _s1(s1) {}
	explicit GF31(const uint32_t n0, const uint32_t n1) : _s0(n0), _s1(n1) {}
	explicit GF31(const uint64_t n0, const uint64_t n1) : _s0(n0), _s1(n1) {}

	const Z31 & s0() const { return _s0; }
	const Z31 & s1() const { return _s1; }

	void set0(const uint32_t n0) { _s0 = Z31(n0); }
	void set1(const uint32_t n1) { _s1 = Z31(n1); }

	bool operator!=(const GF31 & rhs) const { return ((_s0 != rhs._s0) || (_s1 != rhs._s1)); }

	// GF31 conj() const { return GF31(_s0, _s1.neg()); }
	// GF31 muli() const { return GF31(_s1.neg(), _s0); }
	// GF31 half() const { return GF31(_s0.half(), _s1.half()); }

	GF31 operator+(const GF31 & rhs) const { return GF31(_s0 + rhs._s0, _s1 + rhs._s1); }
	GF31 operator-(const GF31 & rhs) const { return GF31(_s0 - rhs._s0, _s1 - rhs._s1); }
	GF31 addconj(const GF31 & rhs) const { return GF31(_s0 + rhs._s0, _s1 - rhs._s1); }
	GF31 subconj(const GF31 & rhs) const { return GF31(_s0 - rhs._s0, _s1 + rhs._s1); }
	GF31 sub_conj(const GF31 & rhs) const { return GF31(_s0 - rhs._s0, rhs._s1 - _s1); }
	GF31 addi(const GF31 & rhs) const { return GF31(_s0 - rhs._s1, _s1 + rhs._s0); }
	GF31 subi(const GF31 & rhs) const { return GF31(_s0 + rhs._s1, _s1 - rhs._s0); }

	GF31 sqr() const { const Z31 t = _s0 * _s1; return GF31(_s0.sqr() - _s1.sqr(), t + t); }
	GF31 mul(const GF31 & rhs) const { return GF31(_s0 * rhs._s0 - _s1 * rhs._s1, _s1 * rhs._s0 + _s0 * rhs._s1); }
	GF31 mulconj(const GF31 & rhs) const { return GF31(_s0 * rhs._s0 + _s1 * rhs._s1, _s1 * rhs._s0 - _s0 * rhs._s1); }

	GF31 lshift(const uint8_t ls0, const uint8_t ls1) const { return GF31(_s0.lshift(ls0), _s1.lshift(ls1)); }
	GF31 rshift(const uint8_t rs0, const uint8_t rs1) const { return GF31(_s0.rshift(rs0), _s1.rshift(rs1)); }

	GF31 pow(const uint64_t e) const
	{
		if (e == 0) return GF31(1u, 0u);
		GF31 r = GF31(1u, 0u), y = *this;
		for (uint64_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r = r.mul(y); y = y.sqr(); }
		return r.mul(y);
	}

	static const GF31 root_nth(const size_t n) { return GF31(Z31(_h_0), Z31(_h_1)).pow(_h_order / n); }
	static uint8_t log2_root_two(const size_t n) { return uint8_t(((uint64_t(1) << 30) / n) % 31); }
};

// Z/{2^61 - 1}Z x Z/{2^31 - 1}Z
class Z61_31
{
private:
	Z61 _n61;
	Z31 _n31;

public:
	Z61_31() {}
	explicit Z61_31(const Z61 & n61, const Z31 & n31) : _n61(n61), _n31(n31) {}

	const Z61 & n61() const { return _n61; }
	const Z31 & n31() const { return _n31; }

	Z61_31 operator+(const Z61_31 & rhs) const { return Z61_31(_n61 + rhs._n61, _n31 + rhs._n31); }
	Z61_31 operator-(const Z61_31 & rhs) const { return Z61_31(_n61 - rhs._n61, _n31 - rhs._n31); }
	Z61_31 operator*(const Z61_31 & rhs) const { return Z61_31(_n61 * rhs._n61, _n31 * rhs._n31); }

	Z61_31 sqr() const { return Z61_31(_n61.sqr(), _n31.sqr()); }
};

struct IBWeight
{
	uint8_t _w61, _w31;

	IBWeight() {}
	IBWeight(const uint8_t w61, const uint8_t w31) : _w61(w61), _w31(w31) {}
	IBWeight operator+(const uint8_t rhs) const { return IBWeight(_w61 + rhs, _w31 + rhs); }
};

// GF((2^61 - 1)^2) x GF((2^31 - 1)^2)
class GF61_31
{
private:
	GF61 _n61;
	GF31 _n31;

public:
	GF61_31() {}
	explicit GF61_31(const uint32_t n) : _n61(GF61(n, 0u)), _n31(GF31(n, 0u)) {}
	explicit GF61_31(const uint64_t n0, const uint64_t n1) : _n61(n0, n1), _n31(n0, n1) {}
	explicit GF61_31(const GF61 & n61, const GF31 & n31) : _n61(n61), _n31(n31) {}
	explicit GF61_31(const Z61_31 & s0, const Z61_31 & s1) : _n61(s0.n61(), s1.n61()), _n31(s0.n31(), s1.n31()) {}

	bool operator!=(const GF61_31 & rhs) const { return ((_n61 != rhs._n61) || (_n31 != rhs._n31)); }

	const Z61_31 s0() const { return Z61_31(_n61.s0(), _n31.s0()); }
	const Z61_31 s1() const { return Z61_31(_n61.s1(), _n31.s1()); }

	// GF61_31 conj() const { return GF61_31(_n61.conj(), _n31.conj()); }
	// GF61_31 muli() const { return GF61_31(_n61.muli(), _n31.muli()); }
	// GF61_31 half() const { return GF61_31(_n61.half(), _n31.half()); }

	GF61_31 operator+(const GF61_31 & rhs) const { return GF61_31(_n61 + rhs._n61, _n31 + rhs._n31); }
	GF61_31 operator-(const GF61_31 & rhs) const { return GF61_31(_n61 - rhs._n61, _n31 - rhs._n31); }
	GF61_31 addconj(const GF61_31 & rhs) const { return GF61_31(_n61.addconj(rhs._n61), _n31.addconj(rhs._n31)); }
	GF61_31 subconj(const GF61_31 & rhs) const { return GF61_31(_n61.subconj(rhs._n61), _n31.subconj(rhs._n31)); }
	GF61_31 sub_conj(const GF61_31 & rhs) const { return GF61_31(_n61.sub_conj(rhs._n61), _n31.sub_conj(rhs._n31)); }
	GF61_31 addi(const GF61_31 & rhs) const { return GF61_31(_n61.addi(rhs._n61), _n31.addi(rhs._n31)); }
	GF61_31 subi(const GF61_31 & rhs) const { return GF61_31(_n61.subi(rhs._n61), _n31.subi(rhs._n31)); }

	GF61_31 sqr() const { return GF61_31(_n61.sqr(), _n31.sqr()); }
	GF61_31 mul(const GF61_31 & rhs) const { return GF61_31(_n61.mul(rhs._n61), _n31.mul(rhs._n31)); }
	GF61_31 mulconj(const GF61_31 & rhs) const { return GF61_31(_n61.mulconj(rhs._n61), _n31.mulconj(rhs._n31)); }

	GF61_31 lshift(const IBWeight ls0, const IBWeight ls1) const { return GF61_31(_n61.lshift(ls0._w61, ls1._w61), _n31.lshift(ls0._w31, ls1._w31)); }
	GF61_31 rshift(const IBWeight rs0, const IBWeight rs1) const { return GF61_31(_n61.rshift(rs0._w61, rs1._w61), _n31.rshift(rs0._w31, rs1._w31)); }

	GF61_31 pow(const uint64_t e) const { return GF61_31(_n61.pow(e), _n31.pow(e)); }

	static const GF61_31 root_nth(const size_t n) { return GF61_31(GF61::root_nth(n), GF31::root_nth(n)); }

	// Chinese remainder theorem
	void garner(__uint128_t & n_0, __uint128_t & n_1) const
	{
		const uint32_t n31_0 = _n31.s0().get(), n31_1 = _n31.s1().get();
		const GF61 n31 = GF61(n31_0, n31_1);
		GF61 u = _n61 - n31; 
		// The inverse of 2^31 - 1 mod 2^61 - 1 is 2^31 + 1
		u = u + u.lshift(31, 31);
		const uint64_t s_0 = u.s0().get(), s_1 = u.s1().get();
		n_0 = n31_0 + (__uint128_t(s_0) << 31) - s_0;
		n_1 = n31_1 + (__uint128_t(s_1) << 31) - s_1;
	}
};

class mersenne
{
private:
	const uint8_t _ln;
	const size_t _n;
	GF61_31 * const _z;
	GF61_31 * const _w;
	IBWeight * const _w_ib;
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
		// the condition is n * (2^{w + 1} - 1)^2 < (2^61 - 1)*(2^31 - 1))
		} while (ln + 2 * (w + 1) >= 61 + 31);

		return ln;
	}

	// Bruun's recursive polynomial factorization
	// See G. Bruun, "z-transform DFT filters and FFT's," in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 26, no. 1, pp. 56-63, February 1978.

	/*void forward2(const size_t m, const size_t s) const
	{
		GF61_31 * const z = _z;
		const GF61_31 * const w = _w;

		for (size_t j = 0; j < s; ++j)
		{
			const GF61_31 wj = w[s + j];
			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 2 * m * j + i;
				const GF61_31 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mul(wj);
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = u0 - u1;
			}
		}
	}*/

	void forward4(const size_t m, const size_t s) const
	{
		GF61_31 * const z = _z;
		const GF61_31 * const w = _w;
		const GF61_31 * const wr4 = &_w[_n / 2];

		for (size_t j = 0; j < s; ++j)
		{
			const GF61_31 w1 = w[s + j], w2 = w[2 * (s + j)], w3 = wr4[s + j];

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 4 * m * j + i;
				const GF61_31 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mul(w2), u2 = z[k + 2 * m].mul(w1), u3 = z[k + 3 * m].mul(w3);
				const GF61_31 v0 = u0 + u2, v1 = u1 + u3, v2 = u0 - u2, v3 = u1 - u3;
				z[k + 0 * m] = v0 + v1; z[k + 1 * m] = v0 - v1;
				z[k + 2 * m] = v2.addi(v3); z[k + 3 * m] = v2.subi(v3);
			}
		}
	}

	// Inverse transform of Bruun's method

	/*void backward2(const size_t m, const size_t s) const
	{
		GF61_31 * const z = _z;
		const GF61_31 * const w = _w;

		for (size_t j = 0; j < s; ++j)
		{
			const GF61_31 wj = w[s + j];
			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 2 * m * j + i;
				const GF61_31 u0 = z[k + 0 * m], u1 = z[k + 1 * m];
				z[k + 0 * m] = u0 + u1;
				z[k + 1 * m] = (u0 - u1).mulconj(wj);
			}
		}
	}*/

	void backward4(const size_t m, const size_t s) const
	{
		GF61_31 * const z = _z;
		const GF61_31 * const w = _w;
		const GF61_31 * const wr4 = &_w[_n / 2];

		for (size_t j = 0; j < s; ++j)
		{
			const GF61_31 w1 = w[s + j], w2 = w[2 * (s + j)], w3 = wr4[s + j];

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 4 * m * j + i;
				const GF61_31 u0 = z[k + 0 * m], u1 = z[k + 1 * m], u2 = z[k + 2 * m], u3 = z[k + 3 * m];
				const GF61_31 v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = u3 - u2;
				z[k + 0 * m] = v0 + v2; z[k + 2 * m] = (v0 - v2).mulconj(w1);
				z[k + 1 * m] = v1.addi(v3).mulconj(w2); z[k + 3 * m] = v1.subi(v3).mulconj(w3);
			}
		}
	}

	// IBDWT: weighted digits
	void weight() const
	{
		GF61_31 * const z = _z;
		const IBWeight * const w_ib = _w_ib;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
			z[k] = z[k].lshift(w_ib[2 * k + 0], w_ib[2 * k + 1]);
		}
	}

	// IBDWT: restore the unweighted digits. sqr and backward transform must be divided by 2n
	void unweight_norm() const
	{
		const uint8_t ln1 = _ln + 1;
		GF61_31 * const z = _z;
		const IBWeight * const w_ib = _w_ib;

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
		GF61_31 * const z = _z;
		const GF61_31 * const w = &_w[_n / 4];

		for (size_t j = 0, n_4 = _n / 4; j < n_4; ++j)
		{
			// const size_t k = 2 * j, kr = bitrev(k, _n / 2), mk = bitrev(_n / 2 - kr, _n / 2);
			const size_t k = 2 * j, mk = (k != 0) ? (size_t(3) << (63 - __builtin_clzll((unsigned long long)k))) - k - 1 : 0;
			const GF61_31 zk = z[k], zmk = z[mk];
			const GF61_31 u0 = zk.addconj(zmk), u1 = zk.subconj(zmk);
			const GF61_31 v0 = u0.sqr() - u1.sqr().mul(w[j]), v1 = u0.mul(u1 + u1);
			z[k] = v0 + v1;
			if (k == 0) z[1] = (z[1] + z[1]).sqr();
			else z[mk] = v0.sub_conj(v1);
		}
	}

	// forward2, sqr, backward2
	void sqr2() const
	{
		const size_t n_4 = _n / 4;
		GF61_31 * const z = _z;
		const GF61_31 * const w = _w;

		for (size_t j = 0; j < n_4; ++j)
		{
			const GF61_31 u0 = z[2 * j + 0], u1 = z[2 * j + 1].mul(w[n_4 + j]);
			z[2 * j + 0] = u0 + u1; z[2 * j + 1] = u0 - u1;
		}

		sqr();

		for (size_t j = 0; j < n_4; ++j)
		{
			const GF61_31 u0 = z[2 * j + 0], u1 = z[2 * j + 1];
			z[2 * j + 0] = u0 + u1; z[2 * j + 1] = (u0 - u1).mulconj(w[n_4 + j]);
		}
	}

	// Add a carry to the number and return the carry of the first digit_width bits
	static constexpr uint64_t digit_adc(const __uint128_t lhs, const uint8_t digit_width, __uint128_t & carry)
	{
		const __uint128_t s = lhs + carry;
		const __uint128_t c = (s < lhs) ? 1 : 0;
		carry = (s >> digit_width) + (c << (128 - digit_width));
		return uint64_t(s) & ((uint64_t(1) << digit_width) - 1);
	}

	// Adjust the digits to the digit representation
	void carry() const
	{
		const size_t n_2 = _n / 2;
		GF61_31 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		__uint128_t c = 0;
		for (size_t k = 0; k < n_2; ++k)
		{
			__uint128_t l0, l1; z[k].garner(l0, l1);
			const uint64_t n0 = digit_adc(l0, digit_width[2 * k + 0], c);
			const uint64_t n1 = digit_adc(l1, digit_width[2 * k + 1], c);
			z[k] = GF61_31(n0, n1);
		}

		while (c != 0)
		{
			for (size_t k = 0; k < n_2; ++k)
			{
				const uint64_t n0 = digit_adc(z[k].s0().n61().get(), digit_width[2 * k + 0], c);
				const uint64_t n1 = digit_adc(z[k].s1().n61().get(), digit_width[2 * k + 1], c);
				z[k] = GF61_31(n0, n1);
				if (c == 0) break;
			}
		}
	}

public:
	mersenne(const uint32_t q) : _ln(transformsize(q)), _n(size_t(1) << _ln),
		_z(new GF61_31[_n / 2]), _w(new GF61_31[_n]), _w_ib(new IBWeight[_n]), _digit_width(new uint8_t[_n])
	{
		const uint8_t ln = _ln;
		const size_t n = _n;

		// radix-2 twiddle factors
		GF61_31 * const w = _w;
		for (size_t s = 1; s <= n / 4; s *= 2)
		{
			const GF61_31 r_s = GF61_31::root_nth(2 * s);
			for (size_t j = 0; j < s; ++j) w[s + j] = r_s.pow(bitrev(j, s));
		}

		// radix-4 twiddle factors
		for (size_t s = 1; s <= n / 4; s *= 2)
		{
			for (size_t j = 0; j < s; ++j) w[n / 2 + s + j] = w[s + j].mul(w[2 * (s + j)]);
		}

		// IBDWT weights: x^q - 1 => x^n - 1
		// See Richard Crandall, Barry Fagin, "Discrete weighted transforms and large-integer arithmetic", Math. Comp. 62 (1994), 305-324.

		// Weights are power of two. Store log_2(weight).
		IBWeight * const w_ib = _w_ib;
		uint8_t * const digit_width = _digit_width;

		const uint8_t lr2_61 = GF61::log2_root_two(n), lr2_31 = GF31::log2_root_two(n);		// n-th root of two

		w_ib[0] = IBWeight(0, 0);

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
					const uint8_t w61 = uint8_t((lr2_61 * (n - r)) % 61);
					const uint8_t w31 = uint8_t((lr2_31 * (n - r)) % 31);
					w_ib[j] = IBWeight(w61, w31);
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
		GF61_31 * const z = _z;

		z[0] = GF61_31(a);
		for (size_t k = 1, n_2 = _n / 2; k < n_2; ++k) z[k] = GF61_31(0u);
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
		GF61_31 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		uint32_t c = a;
		while (c != 0)
		{
			for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
			{
				const uint64_t n0 = digit_sbc(z[k].s0().n61().get(), digit_width[2 * k + 0], c);
				const uint64_t n1 = digit_sbc(z[k].s1().n61().get(), digit_width[2 * k + 1], c);
				z[k] = GF61_31(n0, n1);
				if (c == 0) break;
			}
		}
	}

	bool is_zero() const
	{
		const GF61_31 * const z = _z;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
			if (z[k].s0().n31().get() != 0u) return false;
			if (z[k].s1().n31().get() != 0u) return false;
		}
		return true;
	}

	bool is_Mp() const
	{
		const GF61_31 * const z = _z;
		const uint8_t * const digit_width = _digit_width;

		for (size_t k = 0, n_2 = _n / 2; k < n_2; ++k)
		{
			const uint32_t n0 = z[k].s0().n31().get();
			if (n0 != (uint64_t(1) << digit_width[2 * k + 0]) - 1) return false;
			const uint32_t n1 = z[k].s1().n31().get();
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
