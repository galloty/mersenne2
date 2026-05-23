/*
Copyright 2025, Yves Gallot

mersenne2.cpp is free source code. You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>
#include <array>

inline uint32_t mul_hi(const uint32_t x, const uint32_t y) { return uint32_t((x * uint64_t(y)) >> 32); }
inline uint32_t div3(const uint32_t n) { return mul_hi(n, 1431655766u); }	// = n / 3 if n < 2^31

// Z/M_qZ: the prime field of order M_q = 2^q - 1
// q_inv = 1 (mod q), transform_size | q.
template<int q, uint32_t primroot, uint64_t q_inv, typename uint_t, typename ulong_t>
class ZMq
{
private:
	static const uint_t _p = (uint_t(1) << q) - 1;
	uint_t _n;	// 0 <= n < p

	static uint_t _add(const uint_t a, const uint_t b) { const uint_t t = a + b; return t - ((t >= _p) ? _p : 0); }
	static uint_t _sub(const uint_t a, const uint_t b) { const uint_t t = a - b; return t + ((a < b) ? _p : 0); }
	static uint_t _mul(const uint_t a, const uint_t b) { const ulong_t t = a * ulong_t(b); return _add(uint_t(t) & _p, uint_t(t >> q)); }
	static uint_t _lshift(const uint_t a, const uint8_t s) { const ulong_t t = ulong_t(a) << s; return _add(uint_t(t) & _p, uint_t(t >> q)); }

public:
	ZMq() {}
	explicit ZMq(const uint64_t n) : _n(uint_t((n >= _p) ? n % _p : n)) {}

	uint_t get() const { return _n; }
	static int get_q() { return q; }

	bool operator!=(const ZMq & rhs) const { return (_n != rhs._n); }

	ZMq neg() const { return ZMq((_n == 0) ? 0 : _p - _n); }

	ZMq & operator*=(const ZMq & rhs) { _n = _mul(_n, rhs._n); return *this; }

	ZMq operator+(const ZMq & rhs) const { return ZMq(_add(_n, rhs._n)); }
	ZMq operator-(const ZMq & rhs) const { return ZMq(_sub(_n, rhs._n)); }
	ZMq operator*(const ZMq & rhs) const { return ZMq(_mul(_n, rhs._n)); }

	ZMq sqr() const { return ZMq(_mul(_n, _n)); }

	ZMq lshift(const uint8_t s) const { const uint8_t s_q = s % q; return (s_q != 0) ? ZMq(_lshift(_n, s_q)) : *this; }
	ZMq rshift(const uint8_t s) const { const uint8_t s_q = s % q; return (s_q != 0) ? ZMq(_lshift(_n, q - s_q)) : *this; }

	ZMq pow(const uint64_t e) const
	{
		if (e == 0) return ZMq(1);

		ZMq r = ZMq(1), y = *this;
		for (uint64_t i = e; i != 1; i /= 2)
		{
			if (i % 2 != 0) r = r * y;
			y = y * y;
		}
		r = r * y;

		return r;
	}

	ZMq invert() const { return pow(_p - 2); }

	static const ZMq root_nth(const size_t n) { return ZMq(primroot).pow((_p - 1) / n); }
	static uint8_t log2_root_two(const size_t n) { return uint8_t((q_inv / n) % q); }

	// Cyclic convolution of three numbers
	static void sqr3(ZMq z[3])
	{
		static const ZMq J = root_nth(3);

		// Radix-3
		const ZMq u0 = z[0], u1 = z[1], u2 = z[2];
		const ZMq td = (u1 - u2) * J;
		const ZMq v0 = u0 + u1 + u2, v1 = u0 - u2 + td, v2 = u0 - u1 - td;

		// Squaring
		const ZMq s0 = v0.sqr(), s1 = v1.sqr(), s2 = v2.sqr();

		// Inverse Radix-3
		const ZMq ti = (s1 - s2) * J;
		z[0] = s0 + s1 + s2; z[1] = s0 - s1 - ti; z[2] = s0 - s2 + ti;
	}

	// Weighted convolution of three numbers
	static void weighted_sqr3(ZMq z[3], const ZMq & w, const ZMq & wi)
	{
		// Weighted inputs
		z[1] *= w; z[2] *= w.sqr();

		sqr3(z);

		// Unweighted outputs
		z[1] *= wi; z[2] *= wi.sqr();
	}

	// IBDWT: weighted digits
	static void weight(ZMq * const x, const uint8_t * const w_ib, const size_t n)
	{
		for (size_t k = 0; k < n; ++k) x[k] = x[k].lshift(w_ib[k]);
	}

	// IBDWT: restore the unweighted digits. sqr and backward transform must be divided by 2n or 12n
	static void unweight_norm(ZMq * const x, const uint8_t * const w_ib, const size_t n)
	{
		size_t r = n; uint8_t ln = 1; while (r % 2 == 0) { r /= 2; ++ln; }

		if (r == 1)	// n = 2^m
		{
			for (size_t k = 0; k < n; ++k) x[k] = x[k].rshift(w_ib[k] + ln);
		}
		else		// n = 3 * 2^m
		{
			static const ZMq inv6 = ZMq(6).invert();
			for (size_t k = 0; k < n; ++k) x[k] = x[k].rshift(w_ib[k] + ln) * inv6;
		}
	}
};

// GF(p^2): the field of order p^2, p % 4 = 3
// primroot must be a primitive root of order primroot_order and
//  - the 4th root of unity must be (0, 1),
//  - the 8th root of unity must be 2^{(q - 1)/2} * (1, 1), where p = 2^q - 1.
template<typename Zp, uint64_t primroot_order, uint64_t primroot_0, uint64_t primroot_1>
class GFp2
{
private:
	Zp _s0, _s1;

public:
	GFp2() {}
	explicit GFp2(const Zp & s0, const Zp & s1) : _s0(s0), _s1(s1) {}

	void store(Zp & s0, Zp & s1) const { s0 = _s0; s1 = _s1; }

private:
	GFp2 conj() const { return GFp2(_s0, _s1.neg()); }
	// GFp2 muli() const { return GFp2(_s1.neg(), _s0); }

	GFp2 operator+(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s0, _s1 + rhs._s1); }
	GFp2 operator-(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s0, _s1 - rhs._s1); }
	GFp2 addconj(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s0, _s1 - rhs._s1); }
	GFp2 subconj(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s0, _s1 + rhs._s1); }
	GFp2 sub_conj(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s0, rhs._s1 - _s1); }
	GFp2 addi(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s1, _s1 + rhs._s0); }
	GFp2 subi(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s1, _s1 - rhs._s0); }
	GFp2 addiconj(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s1, _s1 + rhs._s0); }
	GFp2 subiconj(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s1, rhs._s0 - _s1); }
	GFp2 subi_conj(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s1, rhs._s0 - _s1); }

	GFp2 operator*(const Zp & rhs) const { return GFp2(_s0 * rhs, _s1 * rhs); }

	GFp2 sqr() const { const Zp t = _s0 * _s1; return GFp2(_s0.sqr() - _s1.sqr(), t + t); }
	GFp2 mul(const GFp2 & rhs) const { return GFp2(_s0 * rhs._s0 - _s1 * rhs._s1, _s1 * rhs._s0 + _s0 * rhs._s1); }
	GFp2 mulconj(const GFp2 & rhs) const { return GFp2(_s0 * rhs._s0 + _s1 * rhs._s1, _s1 * rhs._s0 - _s0 * rhs._s1); }
	GFp2 mul_neg0(const GFp2 & rhs) const { return GFp2(_s1 * rhs._s1 - _s0 * rhs._s0, _s1 * rhs._s0 + _s0 * rhs._s1); }

	// * sqrt(2)/2 * (1 + i)
	GFp2 mul_R8() const
	{
		const uint8_t log2_sqrt2_2 = uint8_t((Zp::get_q() - 1) / 2);
		return GFp2((_s0 - _s1).lshift(log2_sqrt2_2), (_s1 + _s0).lshift(log2_sqrt2_2));
	}

	// * sqrt(2)/2 * (1 - i)
	GFp2 mul_R8conj() const
	{
		const uint8_t log2_sqrt2_2 = uint8_t((Zp::get_q() - 1) / 2);
		return GFp2((_s0 + _s1).lshift(log2_sqrt2_2), (_s1 - _s0).lshift(log2_sqrt2_2));
	}

	// Radix-2 DIF butterfly
	static void fwd2(Zp * const x0, Zp * const x1, const GFp2 & w, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]);
		const GFp2 u0 = z0 + z1, u1 = (z0 - z1).mul(w);
		u0.store(x0[0], x0[step]); u1.store(x1[0], x1[step]);
	}

	// Radix-2 DIT butterfly
	static void bck2(Zp * const x0, Zp * const x1, const GFp2 & w, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]).mulconj(w);
		const GFp2 u0 = z0 + z1, u1 = z0 - z1;
		u0.store(x0[0], x0[step]); u1.store(x1[0], x1[step]);
	}

	// Radix-4 DIF butterfly
	static void fwd4(Zp * const x0, Zp * const x1, Zp * const x2, Zp * const x3, const GFp2 & w1, const GFp2 & w2, const GFp2 & w3, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]), z2 = GFp2(x2[0], x2[step]), z3 = GFp2(x3[0], x3[step]);
		const GFp2 u0 = z0 + z2, u2 = z0 - z2, u1 = z1 + z3, u3 = z1 - z3;
		const GFp2 v0 = u0 + u1, v1 = (u0 - u1).mul(w2), v2 = u2.addi(u3).mul(w1), v3 = u2.subi(u3).mul(w3);
		v0.store(x0[0], x0[step]); v1.store(x1[0], x1[step]); v2.store(x2[0], x2[step]); v3.store(x3[0], x3[step]);
	}

	// Radix-4 DIT butterfly
	static void bck4(Zp * const x0, Zp * const x1, Zp * const x2, Zp * const x3, const GFp2 & w1, const GFp2 & w2, const GFp2 & w3, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]).mulconj(w2), z2 = GFp2(x2[0], x2[step]).mulconj(w1), z3 = GFp2(x3[0], x3[step]).mulconj(w3);
		const GFp2 u0 = z0 + z1, u1 = z0 - z1, u2 = z2 + z3, u3 = z2 - z3;
		const GFp2 v0 = u0 + u2, v2 = u0 - u2, v1 = u1.subi(u3), v3 = u1.addi(u3);
		v0.store(x0[0], x0[step]); v1.store(x1[0], x1[step]); v2.store(x2[0], x2[step]); v3.store(x3[0], x3[step]);
	}

	// Radix-2 butterfly, last DIF stage or first DIT stage
	static void bty2_0(Zp * const x0, Zp * const x1, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]);
		const GFp2 u0 = z0 + z1, u1 = z0 - z1;
		u0.store(x0[0], x0[step]); u1.store(x1[0], x1[step]);
	}

	// Radix-4 DIF butterfly, last stage
	static void fwd4_0(Zp * const x0, Zp * const x1, Zp * const x2, Zp * const x3, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]), z2 = GFp2(x2[0], x2[step]), z3 = GFp2(x3[0], x3[step]);
		const GFp2 u0 = z0 + z2, u2 = z0 - z2, u1 = z1 + z3, u3 = z1 - z3;
		const GFp2 v0 = u0 + u1, v1 = u0 - u1, v2 = u2.addi(u3), v3 = u2.subi(u3);
		v0.store(x0[0], x0[step]); v1.store(x1[0], x1[step]); v2.store(x2[0], x2[step]); v3.store(x3[0], x3[step]);
	}

	// Radix-4 DIT butterfly, first stage
	static void bck4_0(Zp * const x0, Zp * const x1, Zp * const x2, Zp * const x3, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]), z2 = GFp2(x2[0], x2[step]), z3 = GFp2(x3[0], x3[step]);
		const GFp2 u0 = z0 + z1, u1 = z0 - z1, u2 = z2 + z3, u3 = z2 - z3;
		const GFp2 v0 = u0 + u2, v2 = u0 - u2, v1 = u1.subi(u3), v3 = u1.addi(u3);
		v0.store(x0[0], x0[step]); v1.store(x1[0], x1[step]); v2.store(x2[0], x2[step]); v3.store(x3[0], x3[step]);
	}

	// Radix-8 DIF butterfly, last stage
	static void fwd8_0(Zp * const x0, Zp * const x1, Zp * const x2, Zp * const x3, Zp * const x4, Zp * const x5, Zp * const x6, Zp * const x7, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]), z2 = GFp2(x2[0], x2[step]), z3 = GFp2(x3[0], x3[step]);
		const GFp2 z4 = GFp2(x4[0], x4[step]), z5 = GFp2(x5[0], x5[step]), z6 = GFp2(x6[0], x6[step]), z7 = GFp2(x7[0], x7[step]);

		const GFp2 u0 = z0 + z4, u4 = z0 - z4, u2 = z2 + z6, u6 = z2 - z6;
		const GFp2 u1 = z1 + z5, u5 = z1 - z5, u3 = z3 + z7, u7 = z3 - z7;

		const GFp2 v0 = u0 + u2, v2 = u0 - u2, v4 = u4.addi(u6), v6 = u4.subi(u6);
		const GFp2 v1 = u1 + u3, v3 = u1 - u3, v5 = u5.addi(u7).mul_R8(), v7 = u5.subi(u7).mul_R8();

		const GFp2 t0 = v0 + v1, t1 = v0 - v1, t2 = v2.addi(v3), t3 = v2.subi(v3);
		const GFp2 t4 = v4 + v5, t5 = v4 - v5, t6 = v6.addi(v7), t7 = v6.subi(v7);

		t0.store(x0[0], x0[step]); t1.store(x1[0], x1[step]); t2.store(x2[0], x2[step]); t3.store(x3[0], x3[step]);
		t4.store(x4[0], x4[step]); t5.store(x5[0], x5[step]); t6.store(x6[0], x6[step]); t7.store(x7[0], x7[step]);
	}

	// Radix-8 DIT butterfly, first stage
	static void bck8_0(Zp * const x0, Zp * const x1, Zp * const x2, Zp * const x3, Zp * const x4, Zp * const x5, Zp * const x6, Zp * const x7, const size_t step)
	{
		const GFp2 z0 = GFp2(x0[0], x0[step]), z1 = GFp2(x1[0], x1[step]), z2 = GFp2(x2[0], x2[step]), z3 = GFp2(x3[0], x3[step]);
		const GFp2 z4 = GFp2(x4[0], x4[step]), z5 = GFp2(x5[0], x5[step]), z6 = GFp2(x6[0], x6[step]), z7 = GFp2(x7[0], x7[step]);

		const GFp2 u0 = z0 + z1, u1 = z0 - z1, u2 = z2 + z3, u3 = z2 - z3;
		const GFp2 u4 = z4 + z5, u5 = (z4 - z5).mul_R8conj(), u6 = z6 + z7, u7 = (z6 - z7).mul_R8conj();

		const GFp2 v0 = u0 + u2, v2 = u0 - u2, v4 = u4 + u6, v6 = u4 - u6;
		const GFp2 v1 = u1.subi(u3), v3 = u1.addi(u3), v5 = u5.subi(u7), v7 = u5.addi(u7);

		const GFp2 t0 = v0 + v4, t4 = v0 - v4, t2 = v2.subi(v6), t6 = v2.addi(v6);
		const GFp2 t1 = v1 + v5, t5 = v1 - v5, t3 = v3.subi(v7), t7 = v3.addi(v7);

		t0.store(x0[0], x0[step]); t1.store(x1[0], x1[step]); t2.store(x2[0], x2[step]); t3.store(x3[0], x3[step]);
		t4.store(x4[0], x4[step]); t5.store(x5[0], x5[step]); t6.store(x6[0], x6[step]); t7.store(x7[0], x7[step]);
	}

	// Radix-2 DIF, one iteration; if L = 3 then each point is a vector of three GFp2
	template<size_t L>
	static void forward2(Zp * const x, const GFp2 * const w, const size_t m, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const size_t i_l = (L == 3) ? div3(uint32_t(i)) : i, l = (L == 3) ? 3 * i_l + i : 2 * i;
				const size_t k = 4 * m * j + l;
				fwd2(&x[k + 0 * m], &x[k + 2 * m], w[i_l * s], L);
			}
		}
	}

	// Radix-2 DIT, one iteration; if L = 3 then each point is a vector of three GFp2
	template<size_t L>
	static void backward2(Zp * const x, const GFp2 * const w, const size_t m, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const size_t i_l = (L == 3) ? div3(uint32_t(i)) : i, l = (L == 3) ? 3 * i_l + i : 2 * i;
				const size_t k = 4 * m * j + l;
				bck2(&x[k + 0 * m], &x[k + 2 * m], w[i_l * s], L);
			}
		}
	}

	// Radix-4 DIF, one iteration; if L = 3 then each point is a vector of three GFp2
	template<size_t L>
	static void forward4(Zp * const x, const GFp2 * const w, const size_t m, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const size_t i_l = (L == 3) ? div3(uint32_t(i)) : i, l = (L == 3) ? 3 * i_l + i : 2 * i;
				const size_t k = 8 * m * j + l;
				fwd4(&x[k + 0 * m], &x[k + 2 * m], &x[k + 4 * m], &x[k + 6 * m], w[1 * i_l * s], w[2 * i_l * s], w[3 * i_l * s], L);
			}
		}
	}

	// Radix-4 DIT, one iteration; if L = 3 then each point is a vector of three GFp2
	template<size_t L>
	static void backward4(Zp * const x, const GFp2 * const w, const size_t m, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const size_t i_l = (L == 3) ? div3(uint32_t(i)) : i, l = (L == 3) ? 3 * i_l + i : 2 * i;
				const size_t k = 8 * m * j + l;
				bck4(&x[k + 0 * m], &x[k + 2 * m], &x[k + 4 * m], &x[k + 6 * m], w[1 * i_l * s], w[2 * i_l * s], w[3 * i_l * s], L);
			}
		}
	}

	// Radix-2, one iteration, last DIF stage or first DIT stage
	template<size_t L>
	static void butterfly2_0(Zp * const x, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < L; ++i)
			{
				bty2_0(&x[(4 * j + 0) * L + i], &x[(4 * j + 2) * L + i], L);
			}
		}
	}

	// Radix-4 DIF, one iteration, last stage
	template<size_t L>
	static void forward4_0(Zp * const x, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < L; ++i)
			{
				fwd4_0(&x[(8 * j + 0) * L + i], &x[(8 * j + 2) * L + i], &x[(8 * j + 4) * L + i], &x[(8 * j + 6) * L + i], L);
			}
		}
	}

	// Radix-4 DIT, one iteration, first stage
	template<size_t L>
	static void backward4_0(Zp * const x, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < L; ++i)
			{
				bck4_0(&x[(8 * j + 0) * L + i], &x[(8 * j + 2) * L + i], &x[(8 * j + 4) * L + i], &x[(8 * j + 6) * L + i], L);
			}
		}
	}

	// Radix-8 DIF, one iteration, last stage
	template<size_t L>
	static void forward8_0(Zp * const x, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < L; ++i)
			{
				fwd8_0(&x[(16 * j + 0) * L + i], &x[(16 * j +  2) * L + i], &x[(16 * j +  4) * L + i], &x[(16 * j +  6) * L + i],
					   &x[(16 * j + 8) * L + i], &x[(16 * j + 10) * L + i], &x[(16 * j + 12) * L + i], &x[(16 * j + 14) * L + i], L);
			}
		}
	}

	// Radix-8 DIT, one iteration, first stage
	template<size_t L>
	static void backward8_0(Zp * const x, const size_t s)
	{
		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < L; ++i)
			{
				bck8_0(&x[(16 * j + 0) * L + i], &x[(16 * j +  2) * L + i], &x[(16 * j +  4) * L + i], &x[(16 * j +  6) * L + i],
					   &x[(16 * j + 8) * L + i], &x[(16 * j + 10) * L + i], &x[(16 * j + 12) * L + i], &x[(16 * j + 14) * L + i], L);
			}
		}
	}

	// Weighted convolution of three numbers
	static void weighted_sqr3(GFp2 z[3], const GFp2 & w, const GFp2 & wi)
	{
		static const Zp J = Zp::root_nth(3);

		// Weighted inputs
		const GFp2 u0 = z[0], u1 = z[1].mul(w), u2 = z[2].mul(w.sqr());

		// Radix-3
		const GFp2 td = (u1 - u2) * J;
		const GFp2 v0 = u0 + u1 + u2, v1 = u0 - u2 + td, v2 = u0 - u1 - td;

		// Squaring
		const GFp2 s0 = v0.sqr(), s1 = v1.sqr(), s2 = v2.sqr();

		// Inverse Radix-3
		const GFp2 ti = (s1 - s2) * J;
		const GFp2 t0 = s0 + s1 + s2, t1 = s0 - s1 - ti, t2 = s0 - s2 + ti;

		// Unweighted outputs
		z[0] = t0; z[1] = t1.mul(wi); z[2] = t2.mul(wi.sqr());
	}

	// Input of the transform is 'real' (if alpha is the symbolic square root of p - 1 and the elements of GF(p^2) are a + b*alpha then b = 0).
	// A length-n/2 transform is computed onto z(k) = x(2*k + 0) + x(2*k + 1)*alpha.
	// See Henrik V. Sorensen, Douglas L. Jones, Michael T. Heideman, C. Sidney Burrus, "Real-Valued Fast Fourier Transform Algorithms",
	// in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 35, no. 6, pp. 849-863, June 1987.
	// Output is recombined to produce the half-length transform (full length is not needed because of Hermitian symmetry).
	static void sqr2(Zp * const x, const GFp2 * const w, const size_t n)
	{
		// First point: z[0] = z[-0] and w[0] = 1
		const Zp u0 = x[0] + x[0], u1 = x[1] + x[1];
		x[0] = u0.sqr() + u1.sqr(); x[1] = u0 * (u1 + u1);
		// Middle point: z[n/4] = z[-n/4] and w[n/4] = -1
		const Zp v0 = x[2] + x[2], v1 = x[3] + x[3];
		x[2] = v0.sqr() - v1.sqr(); x[3] = v0 * (v1 + v1);

		for (size_t j = 1; j < n / 4; ++j)
		{
			// const size_t k = 2 * j, kr = bitrev(k, _n / 2), mk = bitrev(_n / 2 - kr, _n / 2);
			const size_t k = 2 * j, mk = (size_t(3) << (63 - __builtin_clzll((unsigned long long)k))) - k - 1;
			const GFp2 zk = GFp2(x[2 * k + 0], x[2 * k + 1]), zmk = GFp2(x[2 * mk + 0], x[2 * mk + 1]);
			const GFp2 u0 = zk.addconj(zmk), u1 = zk.subconj(zmk);
			const GFp2 s0 = u0.sqr() - u1.sqr().mul(w[j]), s1 = u0.mul(u1 + u1);
			const GFp2 v0 = s0 + s1, v1 = s0.sub_conj(s1);
			v0.store(x[2 * k + 0], x[2 * k + 1]); v1.store(x[2 * mk + 0], x[2 * mk + 1]);
		}
	}

	// A "Four Step" NTT algorithm, where n_1 = 2^m and n_2 = 3.
	// See Bailey, D.H. FFTs in external or hierarchical memory. J Supercomput 4, 23–35 (1990).
	// Compute the last stage of the length-n/2 transform (combine z[k] and z[-k] with a {n/3}th root).
	// Multiply by the four-step twiddle factors, calculate the cyclic convolution of length three and apply the reverse operations.
	static void sqr3(Zp * const x, const GFp2 * const w, const size_t n)
	{
		static const Zp R6 = Zp::root_nth(6), R6i = R6.invert();
		static const GFp2 R12 = GFp2::root_nth(12), R12i = R12.invert();

		// First point: z[0] = z[-0] and w[0] = 1
		Zp s0[3], s1[3];
		for (size_t l = 0; l < 3; ++l)
		{
			const Zp u0 = x[0 + l] + x[0 + l], u1 = x[3 + l] + x[3 + l];
			s0[l] = u0 + u1; s1[l] = u0 - u1;
		}
		Zp::sqr3(s0);
		Zp::weighted_sqr3(s1, R6, R6i);
		for (size_t l = 0; l < 3; ++l)
		{
			x[0 + l] = s0[l] + s1[l]; x[3 + l] = s0[l] - s1[l];
		}

		// Middle point: z[n/12] = z[-n/12] and w[n/12] = -1
		GFp2 s2[3];
		for (size_t l = 0; l < 3; ++l) s2[l] = GFp2(x[6 + l] + x[6 + l], x[9 + l] + x[9 + l]);
		weighted_sqr3(s2, R12, R12i);
		for (size_t l = 0; l < 3; ++l) (s2[l] + s2[l]).store(x[6 + l], x[9 + l]);

		for (size_t j = 1; j < n / 12; ++j)
		{
			// const size_t k = 2 * j, kr = bitrev(k, _n / 6), mk = bitrev(_n / 6 - kr, _n / 6);
			const size_t k = 2 * j, mk = (size_t(3) << (63 - __builtin_clzll((unsigned long long)k))) - k - 1;

			GFp2 s0[3], s1[3];
			for (size_t l = 0; l < 3; ++l)
			{
				const GFp2 zk = GFp2(x[6 * k + 0 + l], x[6 * k + 3 + l]), zmk = GFp2(x[6 * mk + 0 + l], x[6 * mk + 3 + l]);
				const GFp2 u0 = zk.addconj(zmk), u1 = zk.subconj(zmk).mul_neg0(w[5 * j + 0]);
				s0[l] = u0.addiconj(u1); s1[l] = u0.subiconj(u1);
			}

			weighted_sqr3(s0, w[5 * j + 1], w[5 * j + 2]);
			weighted_sqr3(s1, w[5 * j + 3], w[5 * j + 4]);

			for (size_t l = 0; l < 3; ++l)
			{
				const GFp2 v0 = s0[l].addconj(s1[l]), v1 = s0[l].subconj(s1[l]).mulconj(w[5 * j + 0]);
				const GFp2 zk = v0.addi(v1), zmk = v0.subi_conj(v1);
				zk.store(x[6 * k + 0 + l], x[6 * k + 3 + l]); zmk.store(x[6 * mk + 0 + l], x[6 * mk + 3 + l]);
			}
		}
	}

public:
	GFp2 pow(const uint64_t e) const
	{
		if (e == 0) return GFp2(Zp(1), Zp(0));
		GFp2 r = GFp2(Zp(1), Zp(0)), y = *this;
		for (uint64_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r = r.mul(y); y = y.sqr(); }
		return r.mul(y);
	}

	GFp2 invert() const
	{
		// 1 / (s0 + i * s1) = (s0 - i * s1) / (s0^2 + s1^2)
		return conj() * (_s0.sqr() + _s1.sqr()).invert();
	}

	static const GFp2 root_nth(const size_t n) { return GFp2(Zp(primroot_0), Zp(primroot_1)).pow(primroot_order / n); }

	static void square_2(Zp * const x, const GFp2 * const w, const size_t n)
	{
		// for (size_t m = n / 2, s = 1; m >= 2; m /= 2, s *= 2) forward2<1>(x, w, m / 2, s);
		size_t m = n / 2, s = 1; for (; m > 8; m /= 4, s *= 4) forward4<1>(x, w, m / 4, s);
		if (m == 8) forward8_0<1>(x, n / 16);
		else if (m == 4) forward4_0<1>(x, n / 8);
		else if (m == 2) butterfly2_0<1>(x, n / 4);

		sqr2(x, &w[n / 2], n);

		// for (size_t m = 2, s = n / 4; s >= 1; m *= 2, s /= 2) backward2<1>(x, w, m / 2, s);
		if (m == 8) backward8_0<1>(x, n / 16);
		else if (m == 4) backward4_0<1>(x, n / 8);
		else if (m == 2) butterfly2_0<1>(x, n / 4);
		m *= 4; s /= 4; for (; s >= 1; m *= 4, s /= 4) backward4<1>(x, w, m / 4, s);

	}

	static void square_3(Zp * const x, const GFp2 * const w, const size_t n)
	{
		// for (size_t m = n / 2, s = 1; m >= 6; m /= 2, s *= 2) forward2<3>(x, w, m / 2, s);
		size_t m = n / 2, s = 1; for (; m > 24; m /= 4, s *= 4) forward4<3>(x, w, m / 4, s);
		if (m == 24) forward8_0<3>(x, n / 48);
		else if (m == 12) forward4_0<3>(x, n / 24);
		else if (m == 6) butterfly2_0<3>(x, n / 12);

		sqr3(x, &w[n / 6], n);

		// for (size_t m = 6, s = n / 12; s >= 1; m *= 2, s /= 2) backward2<3>(x, w, m / 2, s);
		if (m == 24) backward8_0<3>(x, n / 48);
		else if (m == 12) backward4_0<3>(x, n / 24);
		else if (m == 6) butterfly2_0<3>(x, n / 12);
		m *= 4; s /= 4; for (; s >= 1; m *= 4, s /= 4) backward4<3>(x, w, m / 4, s);
	}
};

using Z61 = ZMq<61, 37, uint64_t(3) << 54, uint64_t, __uint128_t>;
using Z31 = ZMq<31, 7, uint64_t(33) << 54, uint32_t, uint64_t>;

using GF61 = GFp2<Z61, uint64_t(9) << 30, 498212544045757ull, 3337356182474291ull>;
using GF31 = GFp2<Z31, uint64_t(9) << 30, 690004u, 2724878u>;

class mersenne
{
private:
	const size_t _n;
	Z61 * const _x61;
	Z31 * const _x31;
	GF61 * const _w61;
	GF31 * const _w31;
	uint8_t * const _w_ib61;
	uint8_t * const _w_ib31;
	uint8_t * const _digit_width;

private:
	// Make sure the transform is long enough so that each digit cannot overflow after the convolution
	static constexpr size_t transform_size(const uint32_t q)
	{
		uint8_t log2_n = 1, log2_n3 = 1; uint32_t w = 0;
		do
		{
			++log2_n;
			// digit-width is w or w + 1
			w = q >> log2_n;
		// the condition is n * (2^{w + 1} - 1)^2 < (2^61 - 1)*(2^31 - 1))
		// 2 * (w + 1) + log2(n) <= 91 < 91.99999999932
		} while (2 * (w + 1) + log2_n > 91);

		do
		{
			++log2_n3;
			w = q / (3u << log2_n3);
		// The condition is 3 * n * (2^{w + 1} - 1)^2 < (2^61 - 1)*(2^31 - 1))
		// 2 * (w + 1) + log2(n) <= 90 < 91.99999999932 - 1.5849625
		} while (2 * (w + 1) + log2_n3 > 90);

		return std::min(size_t(1) << log2_n, size_t(3) << log2_n3);
	}

	// Chinese remainder theorem
	static __uint128_t garner(const Z61 & n61, const Z31 & n31)
	{
		const Z61 u = n61 - Z61(n31.get()); 
		// The inverse of 2^31 - 1 mod 2^61 - 1 is 2^31 + 1
		const Z61 v = u + u.lshift(31);
		const uint64_t s = v.get();
		return n31.get() + (__uint128_t(s) << 31) - s;
	}

	// Add a carry to the number and return the carry of the first digit_width bits
	static constexpr uint64_t digit_adc(const __uint128_t lhs, const uint8_t digit_width, __uint128_t & carry)
	{
		const __uint128_t s = lhs + carry;
		const __uint128_t c = (s < lhs) ? 1 : 0;
		carry = (s >> digit_width) + (c << (128 - digit_width));
		return uint64_t(s) & ((uint64_t(1) << digit_width) - 1);
	}

	// Subtract a carry to the number and return the carry if borrowing
	static constexpr uint64_t digit_sbc(const uint64_t lhs, const uint8_t digit_width, uint32_t & carry)
	{
		const bool borrow = (lhs < carry);
		const uint64_t r = lhs - carry + (borrow ? (uint64_t(1) << digit_width) : 0);
		carry = borrow ? 1 : 0;
		return r;
	}

	// Adjust the digits to the digit representation
	void carry() const
	{
		Z61 * const x61 = _x61;
		Z31 * const x31 = _x31;
		const uint8_t * const digit_width = _digit_width;
		const size_t n = _n;

		__uint128_t c = 0;
		for (size_t k = 0; k < n; ++k)
		{
			const __uint128_t l = garner(x61[k], x31[k]);
			const uint64_t r = digit_adc(l, digit_width[k], c);
			x61[k] = Z61(r); x31[k] = Z31(r);
		}

		while (c != 0)
		{
			for (size_t k = 0; k < n; ++k)
			{
				const uint64_t r = digit_adc(x61[k].get(), digit_width[k], c);
				x61[k] = Z61(r); x31[k] = Z31(r);
				if (c == 0) break;
			}
		}
	}

	// Bit-reversal permutation of index i for a sequence of n items
	static constexpr size_t bitrev(const size_t i, const size_t n)
	{
		size_t r = 0;
		for (size_t k = n, j = i; k > 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
		return r;
	}

	template<typename GF>
	static void init_twiddle_factors(GF * const w, const size_t n)
	{
		// Radix-2
		const size_t n3 = (n % 3 != 0) ? n : n / 3;
		const GF r = GF::root_nth(n3 / 2);
		for (size_t j = 0; j < n3 / 2; ++j) w[j] = r.pow(j);

		// Hermitian product
		GF * const wh = &w[n3 / 2];
		if (n % 3 != 0)	// n = 2^m
		{
			const GF r_n_2 = GF::root_nth(n / 2);
			for (size_t j = 0; j < n / 4; ++j) wh[j] = r_n_2.pow(bitrev(j, n / 4));
		}
		else			// n = 3 * 2^m
		{
			const GF r_n = GF::root_nth(n), r_n_3 = GF::root_nth(n / 3);
			for (size_t j = 0; j < n / 12; ++j)
			{
				const size_t jr = bitrev(j, n / 12);
				wh[5 * j + 0] = r_n_3.pow(bitrev(j, n / 12));	// Four step factor
				const GF rj = r_n.pow(jr);
				wh[5 * j + 1] = rj; wh[5 * j + 2] = rj.invert();
				const GF rjp = r_n.pow(n / 6 - jr);
				wh[5 * j + 3] = rjp; wh[5 * j + 4] = rjp.invert();
			}
		}
	}

public:
	mersenne(const uint32_t q) : _n(transform_size(q)),
		_x61(new Z61[_n]), _x31(new Z31[_n]), _w61(new GF61[_n]), _w31(new GF31[_n]),
		_w_ib61(new uint8_t[_n]), _w_ib31(new uint8_t[_n]), _digit_width(new uint8_t[_n])
	{
		const size_t n = _n;

		init_twiddle_factors<GF61>(_w61, n);
		init_twiddle_factors<GF31>(_w31, n);

		// IBDWT weights: x^q - 1 => x^n - 1
		// See Richard Crandall, Barry Fagin, "Discrete weighted transforms and large-integer arithmetic", Math. Comp. 62 (1994), 305-324.

		// Weights are power of two: store log_2(weight)
		uint8_t * const w_ib61 = _w_ib61;
		uint8_t * const w_ib31 = _w_ib31;
		uint8_t * const digit_width = _digit_width;

		const uint8_t lr2_61 = Z61::log2_root_two(n), lr2_31 = Z31::log2_root_two(n);		// n-th root of two

		w_ib61[0] = 0; w_ib31[0] = 0;

		const uint64_t q_n = uint64_t(q / n);

		uint64_t ceil_qjm1_n = 0;
		for (size_t j = 1; j <= n; ++j)
		{
			const uint64_t qj = q * uint64_t(j);
			// ceil(a / b) = floor((a - 1) / b) + 1
			const uint64_t ceil_qj_n = (j == n) ? q : (qj - 1) / n + 1;

			// Bit position for digit[i] is ceil(qj / n)
			const uint64_t c = ceil_qj_n - ceil_qjm1_n;
			if ((c != q_n) && (c != q_n + 1u)) throw;
			digit_width[j - 1] = uint8_t(c);

			if (j == n) break;

			// Weight is 2^[ceil(qj/n) - qj/n]
			// e = (ceil(qj / n).n - qj) / n
			// qj = k.n => e = 0
			// qj = k.n + r, r > 0 => ((k + 1).n - (k.n + r)) / n = (n - r) / n
			const uint64_t r = qj % n;
			w_ib61[j] = uint8_t((r != 0) ? (lr2_61 * (n - r)) % 61 : 0);
			w_ib31[j] = uint8_t((r != 0) ? (lr2_31 * (n - r)) % 31 : 0);

			ceil_qjm1_n = ceil_qj_n;
		}
	}

	virtual ~mersenne()
	{
		delete[] _x61;
		delete[] _x31;
		delete[] _w61;
		delete[] _w31;
		delete[] _w_ib61;
		delete[] _w_ib31;
		delete[] _digit_width;
	}

	size_t get_length() const { return _n; }

	void init(const uint32_t a) const
	{
		Z61 * const x61 = _x61;
		Z31 * const x31 = _x31;
		const size_t n = _n;

		x61[0] = Z61(a); for (size_t k = 1; k < n; ++k) x61[k] = Z61(0);
		x31[0] = Z31(a); for (size_t k = 1; k < n; ++k) x31[k] = Z31(0);
	}

	void square() const
	{
		Z61 * const x61 = _x61;
		Z31 * const x31 = _x31;
		const GF61 * const w61 = _w61;
		const GF31 * const w31 = _w31;
		const uint8_t * const w_ib61 = _w_ib61;
		const uint8_t * const w_ib31 = _w_ib31;
		const size_t n = _n;

		// Weighted convolutions
		if (_n % 3 != 0)	// n = 2^m
		{
			Z61::weight(x61, w_ib61, n);
			GF61::square_2(x61, w61, n);
			Z61::unweight_norm(x61, w_ib61, n);

			Z31::weight(x31, w_ib31, n);
			GF31::square_2(x31, w31, n);
			Z31::unweight_norm(x31, w_ib31, n);
		}
		else				// n = 3 * 2^m
		{
			Z61::weight(x61, w_ib61, n);
			GF61::square_3(x61, w61, n);
			Z61::unweight_norm(x61, w_ib61, n);

			Z31::weight(x31, w_ib31, n);
			GF31::square_3(x31, w31, n);
			Z31::unweight_norm(x31, w_ib31, n);
		}

		// Carry propagation
		carry();
	}

	void sub(const uint32_t a) const
	{
		Z61 * const x61 = _x61;
		Z31 * const x31 = _x31;
		const uint8_t * const digit_width = _digit_width;
		const size_t n = _n;

		uint32_t c = a;
		while (c != 0)
		{
			for (size_t k = 0; k < n; ++k)
			{
				const uint64_t r = digit_sbc(x61[k].get(), digit_width[k], c);
				x61[k] = Z61(r); x31[k] = Z31(r);
				if (c == 0) break;
			}
		}
	}

	bool is_zero() const
	{
		const Z61 * const x61 = _x61;
		const size_t n = _n;

		for (size_t k = 0; k < n; ++k)
		{
			if (x61[k].get() != 0) return false;
		}
		return true;
	}

	bool is_Mp() const
	{
		const Z61 * const x61 = _x61;
		const uint8_t * const digit_width = _digit_width;
		const size_t n = _n;

		for (size_t k = 0; k < n; ++k)
		{
			if (x61[k].get() != (uint64_t(1) << digit_width[k]) - 1) return false;
		}

		return true;
	}
};

#define	CHECK_MERSENNE_PRIMES	true

int main()
{
#ifdef CHECK_MERSENNE_PRIMES
	std::array<uint32_t, 51> mexp = { 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941,
		11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593,
		13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933, 136279841 };
#endif

#ifdef CHECK_MERSENNE_PRIMES
	for (uint32_t p : mexp)
#else
	for (uint32_t p = 3; p <= 4294967291; p += 2)
#endif
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
		if (m.is_zero() || m.is_Mp()) std::cout << p << ", " << m.get_length() << std::endl;
#ifdef CHECK_MERSENNE_PRIMES
		else std::cout << p << " failed!" << std::endl;
#endif
	}

	return EXIT_SUCCESS;
}
