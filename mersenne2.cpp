/*
Copyright 2025, Yves Gallot

mersenne2.cpp is free source code. You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <iostream>
#include <cstdint>
#include <array>

// TODO move to mersenne
// bit-reversal permutation of index i for a sequence of n items
static constexpr size_t bitrev(const size_t i, const size_t n)
{
	size_t r = 0;
	for (size_t k = n, j = i; k > 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
	return r;
}

// Z/M_qZ: the prime field of order M_q = 2^q - 1
// q_inv = 1 (mod q), transform_size | q
template<int q, uint32_t primroot, uint64_t q_inv, typename uint_t, typename ulong_t>
class ZMq
{
private:
	static const uint_t _p = (uint_t(1) << q) - 1;
	uint_t _n;	// 0 <= n < p

	static uint_t _add(const uint_t a, const uint_t b)
	{
		const uint_t t = a + b;
		return t - ((t >= _p) ? _p : 0);
	}

	static uint_t _sub(const uint_t a, const uint_t b)
	{
		const uint_t t = a - b;
		return t + ((a < b) ? _p : 0);
	}

	static uint_t _mul(const uint_t a, const uint_t b)
	{
		const ulong_t t = a * ulong_t(b);
		return _add(uint_t(t) & _p, uint_t(t >> q));
	}

	static uint_t _lshift(const uint_t a, const uint8_t s)
	{
		const ulong_t t = ulong_t(a) << s;
		return _add(uint_t(t) & _p, uint_t(t >> q));
	}

public:
	ZMq() {}
	explicit ZMq(const uint64_t n) : _n(uint_t((n >= _p) ? n % _p : n)) {}

	uint_t get() const { return _n; }
	static uint_t get_p() { return _p; }

	bool operator!=(const ZMq & rhs) const { return (_n != rhs._n); }

	ZMq neg() const { return ZMq((_n == 0) ? 0 : _p - _n); }
	ZMq half() const { return ZMq(((_n % 2 == 0) ? _n : (_n + _p)) / 2); }

	ZMq operator+(const ZMq & rhs) const { return ZMq(_add(_n, rhs._n)); }
	ZMq operator-(const ZMq & rhs) const { return ZMq(_sub(_n, rhs._n)); }
	ZMq operator*(const ZMq & rhs) const { return ZMq(_mul(_n, rhs._n)); }

	ZMq sqr() const { return ZMq(_mul(_n, _n)); }

	ZMq lshift(const uint8_t s) const { const uint8_t s_q = s % q; return (s_q != 0) ? ZMq(_lshift(_n, s_q)) : *this; }
	ZMq rshift(const uint8_t s) const { const uint8_t s_q = s % q; return (s_q != 0) ? ZMq(_lshift(_n, q - s_q)) : *this; }

	ZMq pow(const uint64_t e) const
	{
		if (e == 0) return ZMq(1u);

		ZMq r = ZMq(1u), y = *this;
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

	// IBDWT: weighted digits
	static void weight(ZMq * const x, const uint8_t * const w_ib, const size_t n)
	{
		for (size_t k = 0; k < n; ++k) x[k] = x[k].lshift(w_ib[k]);
	}

	// IBDWT: restore the unweighted digits. sqr and backward transform must be divided by 2n or 12n
	static void unweight_norm(ZMq * const x, const uint8_t * const w_ib, const size_t n)
	{
		size_t r = n; uint8_t ln = 1; while (r % 2 == 0) { r /= 2; ++ln; }

		if (r == 1)
		{
			for (size_t k = 0; k < n; ++k) x[k] = x[k].rshift(w_ib[k] + ln);
		}
		else	// n = 3 * 2^m
		{
			static const ZMq inv6 = ZMq(6).invert();
			for (size_t k = 0; k < n; ++k) x[k] = x[k].rshift(w_ib[k] + ln) * inv6;
		}
	}
};

// GF(p^2): the prime field of order p^2, p % 4 = 3
// primroot must be a primitive root of order primroot_order which is a root of (0, 1).
template<typename Zp, uint64_t primroot_order, uint64_t primroot_0, uint64_t primroot_1>
class GFp2
{
private:
	Zp _s0, _s1;

public:
	GFp2() {}
	explicit GFp2(const Zp & s0, const Zp & s1) : _s0(s0), _s1(s1) {}
	explicit GFp2(const uint64_t n0, const uint64_t n1) : _s0(n0), _s1(n1) {}

	void store(Zp & s0, Zp & s1) const { s0 = _s0; s1 = _s1; }

	const Zp & s0() const { return _s0; }
	const Zp & s1() const { return _s1; }
	Zp & s0() { return _s0; }
	Zp & s1() { return _s1; }

	void set0(const uint64_t n0) { _s0 = Zp(n0); }
	void set1(const uint64_t n1) { _s1 = Zp(n1); }

	bool operator!=(const GFp2 & rhs) const { return ((_s0 != rhs._s0) || (_s1 != rhs._s1)); }

	GFp2 conj() const { return GFp2(_s0, _s1.neg()); }
	GFp2 muli() const { return GFp2(_s1.neg(), _s0); }
	GFp2 half() const { return GFp2(_s0.half(), _s1.half()); }

	GFp2 operator+(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s0, _s1 + rhs._s1); }
	GFp2 operator-(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s0, _s1 - rhs._s1); }
	GFp2 addconj(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s0, _s1 - rhs._s1); }
	GFp2 subconj(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s0, _s1 + rhs._s1); }
	GFp2 sub_conj(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s0, rhs._s1 - _s1); }
	GFp2 addi(const GFp2 & rhs) const { return GFp2(_s0 - rhs._s1, _s1 + rhs._s0); }
	GFp2 subi(const GFp2 & rhs) const { return GFp2(_s0 + rhs._s1, _s1 - rhs._s0); }

	GFp2 operator*(const Zp & rhs) const { return GFp2(_s0 * rhs, _s1 * rhs); }

	GFp2 sqr() const { const Zp t = _s0 * _s1; return GFp2(_s0.sqr() - _s1.sqr(), t + t); }
	GFp2 mul(const GFp2 & rhs) const { return GFp2(_s0 * rhs._s0 - _s1 * rhs._s1, _s1 * rhs._s0 + _s0 * rhs._s1); }
	GFp2 mulconj(const GFp2 & rhs) const { return GFp2(_s0 * rhs._s0 + _s1 * rhs._s1, _s1 * rhs._s0 - _s0 * rhs._s1); }

	GFp2 pow(const __uint128_t e) const
	{
		if (e == 0) return GFp2(1u, 0u);
		GFp2 r = GFp2(1u, 0u), y = *this;
		for (__uint128_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r = r.mul(y); y = y.sqr(); }
		return r.mul(y);
	}

	GFp2 invert() const { const __uint128_t t = Zp::get_p(); return pow(t * t - 2); }

	static const GFp2 root_nth(const size_t n) { return GFp2(Zp(primroot_0), Zp(primroot_1)).pow(primroot_order / n); }

	// Radix-2 DIF
	template<size_t L>
	static void forward2(Zp * const x, const size_t m, const size_t s)
	{
		const GFp2 r_m = root_nth(2 * m);

		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = L * (2 * m * j + i);
				const GFp2 r_mi = r_m.pow(i);

				for (size_t l = 0; l < L; ++l)
				{
					const GFp2 u0 = GFp2(x[2 * (k + 0 * m) + 0 + l], x[2 * (k + 0 * m) + L + l]);
					const GFp2 u1 = GFp2(x[2 * (k + L * m) + 0 + l], x[2 * (k + L * m) + L + l]);
					const GFp2 v0 = u0 + u1, v1 = (u0 - u1).mul(r_mi);
					v0.store(x[2 * (k + 0 * m) + 0 + l], x[2 * (k + 0 * m) + L + l]);
					v1.store(x[2 * (k + L * m) + 0 + l], x[2 * (k + L * m) + L + l]);
				}
			}
		}
	}

	// Radix-2 DIT
	template<size_t L>
	static void backward2(Zp * const x, const size_t m, const size_t s)
	{
		const GFp2 r_m = root_nth(2 * m);

		for (size_t j = 0; j < s; ++j)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = L * (2 * m * j + i);
				const GFp2 r_mi = r_m.pow(i);

				for (size_t l = 0; l < L; ++l)
				{
					const GFp2 u0 = GFp2(x[2 * (k + 0 * m) + 0 + l], x[2 * (k + 0 * m) + L + l]);
					const GFp2 u1 = GFp2(x[2 * (k + L * m) + 0 + l], x[2 * (k + L * m) + L + l]).mulconj(r_mi);
					const GFp2 v0 = u0 + u1, v1 = u0 - u1;
					v0.store(x[2 * (k + 0 * m) + 0 + l], x[2 * (k + 0 * m) + L + l]);
					v1.store(x[2 * (k + L * m) + 0 + l], x[2 * (k + L * m) + L + l]);
				}
			}
		}
	}

	static void weighted_sqr3(GFp2 z[3], const GFp2 w, GFp2 wi)
	{
		const GFp2 w2 = w.sqr(), wi2 = wi.sqr();
		static const Zp J = Zp::root_nth(3);

		// weighted inputs
		const GFp2 u0 = z[0], u1 = z[1].mul(w), u2 = z[2].mul(w2);

		// Radix-3
		const GFp2 td = (u1 - u2) * J;
		const GFp2 v0 = u0 + u1 + u2, v1 = u0 - u2 + td, v2 = u0 - u1 - td;

		// square
		const GFp2 s0 = v0.sqr(), s1 = v1.sqr(), s2 = v2.sqr();

		// Inverse Radix-3
		const GFp2 ti = (s1 - s2) * J;
		const GFp2 t0 = s0 + s1 + s2, t1 = s0 - s1 - ti, t2 = s0 - s2 + ti;

		// unweighted outputs
		z[0] = t0; z[1] = t1.mul(wi); z[2] = t2.mul(wi2);
	}


	// Input of the transform is 'real' (if alpha is the symbolic square root of p - 1 and the elements of GF(p^2) are a + b*alpha then b = 0).
	// A length-n/2 transform is computed onto z(k) = x(2*k + 0) + x(2*k + 1)*alpha.
	// See Henrik V. Sorensen, Douglas L. Jones, Michael T. Heideman, C. Sidney Burrus, "Real-Valued Fast Fourier Transform Algorithms",
	// in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 35, no. 6, pp. 849-863, June 1987.
	// Output is recombined to produce the half-length transform (full length is not needed because of Hermitian symmetry).
	static void sqr(Zp * const x, const GFp2 * const w, const size_t n)
	{
		for (size_t j = 0; j < n / 4; ++j)
		{
			// const size_t k = 2 * j, kr = bitrev(k, _n / 2), mk = bitrev(_n / 2 - kr, _n / 2);
			const size_t k = 2 * j, mk = (k != 0) ? (size_t(3) << (63 - __builtin_clzll((unsigned long long)k))) - k - 1 : 0;
			const GFp2 zk = GFp2(x[2 * k + 0], x[2 * k + 1]), zmk = GFp2(x[2 * mk + 0], x[2 * mk + 1]);
			const GFp2 u0 = zk.addconj(zmk), u1 = zk.subconj(zmk);
			const GFp2 s0 = u0.sqr() - u1.sqr().mul(w[j]), s1 = u0.mul(u1 + u1);
			const GFp2 v0 = s0 + s1, v1 = s0.sub_conj(s1);
			v0.store(x[2 * k + 0], x[2 * k + 1]);
			if (k == 0) { GFp2 z1 = GFp2(x[2], x[3]); z1 = (z1 + z1).sqr(); z1.store(x[2], x[3]); }
			else v1.store(x[2 * mk + 0], x[2 * mk + 1]);
		}
	}

	// A "Four Step" NTT algorithm, where n_1 = 2^m and n_2 = 3.
	// See Bailey, D.H. FFTs in external or hierarchical memory. J Supercomput 4, 23–35 (1990).
	// Compute the last stage of the length-n/2 transform (combine z[k] and z[-k] with a {n/3}th root).
	// Multiply by the four-step twiddle factors. Calculate the radix-3, square and apply the reverse operations.
	static void sqr3(Zp * const x, const size_t n)
	{
		const GFp2 r = root_nth(n / 3);
		const GFp2 w_n = root_nth(n);

		for (size_t j = 0; j < n / 12; ++j)
		{
			// const size_t k = 2 * j, kr = bitrev(k, _n / 6), mk = bitrev(_n / 6 - kr, _n / 6);
			const size_t k = 2 * j, mk = (k != 0) ? (size_t(3) << (63 - __builtin_clzll((unsigned long long)k))) - k - 1 : 0;
			const GFp2 rk = r.pow(bitrev(j, n / 12));
			GFp2 s0[3], s1[3];
			for (size_t l = 0; l < 3; ++l)
			{
				const GFp2 zk = GFp2(x[6 * k + 0 + l], x[6 * k + 3 + l]), zmk = GFp2(x[6 * mk + 0 + l], x[6 * mk + 3 + l]);
				const GFp2 u0 = zk.addconj(zmk), u1 = zk.subconj(zmk).mul(rk).muli();
				s0[l] = u0 - u1; s1[l] = (k == 0) ? u0 + u1 : (u0 + u1).conj();
			}

			const GFp2 w_nj = w_n.pow(bitrev(j, n / 12));
			weighted_sqr3(s0, w_nj, w_nj.invert());
			const GFp2 w_nj_6 = w_n.pow(n / 6 - bitrev(j, n / 12));
			weighted_sqr3(s1, w_nj_6, w_nj_6.invert());

			for (size_t l = 0; l < 3; ++l)
			{
				const GFp2 v0 = s0[l] + s1[l].conj(), v1 = (s0[l] - s1[l].conj()).mul(rk.conj()).muli();
				const GFp2 zk = v0 + v1;
				zk.store(x[6* k + 0 + l], x[6 * k + 3 + l]);
				if (k != 0)
				{
					const GFp2 zmk = (v0 - v1).conj();
					zmk.store(x[6 * mk + 0 + l], x[6 * mk + 3 + l]);
				}
			}

			if (k == 0)
			{
				GFp2 z1[3];
				for (size_t l = 0; l < 3; ++l) z1[l] = GFp2(x[6 + 0 + l], x[6 + 3 + l]);
				for (size_t l = 0; l < 3; ++l) z1[l] = z1[l] + z1[l];

				const GFp2 w_n12 = w_n.pow(n / 12);
				weighted_sqr3(z1, w_n12, w_n12.invert());

				for (size_t l = 0; l < 3; ++l) z1[l] = z1[l] + z1[l];
				for (size_t l = 0; l < 3; ++l) z1[l].store(x[6 + 0 + l], x[6 + 3 + l]);
			}
		}
	}
};

using Z61 = ZMq<61, 37, uint64_t(3) << 54, uint64_t, __uint128_t>;
using Z31 = ZMq<31, 7, uint64_t(33) << 54, uint32_t, uint64_t>;

using GF61 = GFp2<Z61, uint64_t(3) << 62, 1794883824738941ull, 7671273768663199ull>;
using GF31 = GFp2<Z31, uint64_t(3) << 32, 698590u, 9209u>;

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
	// Make sure the transform is long enough so that each digit cannot overflow after the convolution.
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
		} while (2 * (w + 1) * 2 + log2_n > 91);

		do
		{
			++log2_n3;
			w = q / (3u << log2_n3);
		// The condition is 3 * n * (2^{w + 1} - 1)^2 < (2^61 - 1)*(2^31 - 1))
		// 2 * (w + 1) + log2(n) <= 90 < 91.99999999932 - 1.5849625
		} while ((w + 1) * 2 + log2_n3 > 90);

		// return size_t(1) << log2_n;
		// return size_t(3) << log2_n3;
		return std::min(size_t(1) << log2_n, size_t(3) << log2_n3);
	}

	// Chinese remainder theorem
	static __uint128_t garner(const Z61 & n61, const Z31 & n31)
	{
		Z61 u = n61 - Z61(n31.get()); 
		// The inverse of 2^31 - 1 mod 2^61 - 1 is 2^31 + 1
		u = u + u.lshift(31);
		const uint64_t s = u.get();
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

public:
	mersenne(const uint32_t q) : _n(transform_size(q)),
		_x61(new Z61[_n]), _x31(new Z31[_n]),
		_w61(new GF61[2 * _n]), _w31(new GF31[2 * _n]),
		_w_ib61(new uint8_t[_n]), _w_ib31(new uint8_t[_n]),
		_digit_width(new uint8_t[_n])
	{
		// const uint8_t ln = _ln;
		const size_t n = _n;

		// radix-2 twiddle factors
		GF61 * const w61 = _w61;
		GF31 * const w31 = _w31;
		// for (size_t s = 1; s < n / 4; s *= 2)
		// {
		// 	const GF61_31 r_s = GF61_31::root_nth(2 * s);
		// 	for (size_t j = 0; j < s; ++j) w[s + j] = r_s.pow(bitrev(j, s));
		// }
		if (n % 3 != 0)
		{
			// Hermitian product twiddle factors
			size_t s = n / 4;
			const GF61 r_s61 = GF61::root_nth(2 * s);
			for (size_t j = 0; j < s; ++j) w61[s + j] = r_s61.pow(bitrev(j, s));
			const GF31 r_s31 = GF31::root_nth(2 * s);
			for (size_t j = 0; j < s; ++j) w31[s + j] = r_s31.pow(bitrev(j, s));
		}

		// radix-4 twiddle factors
		// for (size_t s = 1; s <= n / 4; s *= 2)
		// {
		// 	for (size_t j = 0; j < s; ++j) w[n / 2 + s + j] = w[s + j].mul(w[2 * (s + j)]);
		// }

		// IBDWT weights: x^q - 1 => x^n - 1
		// See Richard Crandall, Barry Fagin, "Discrete weighted transforms and large-integer arithmetic", Math. Comp. 62 (1994), 305-324.

		// Weights are power of two. Store log_2(weight).
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

			// bit position for digit[i] is ceil(qj / n)
			const uint64_t c = ceil_qj_n - ceil_qjm1_n;
			if ((c != q_n) && (c != q_n + 1u)) throw;
			digit_width[j - 1] = uint8_t(c);

			if (j == n) break;

			// weight is 2^[ceil(qj/n) - qj/n]
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

		x61[0] = Z61(a); for (size_t k = 1; k < n; ++k) x61[k] = Z61(0u);
		x31[0] = Z31(a); for (size_t k = 1; k < n; ++k) x31[k] = Z31(0u);
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

		if (_n % 3 != 0)
		{
			// weighted convolution, radix-2 transforms

			Z61::weight(x61, w_ib61, n);
			for (size_t m = n / 4, s = 1; m >= 1; m /= 2, s *= 2) GF61::forward2<1>(x61, m, s);
			GF61::sqr(x61, &w61[n / 4], n);
			for (size_t m = 1, s = n / 4; s >= 1; m *= 2, s /= 2) GF61::backward2<1>(x61, m, s);
			Z61::unweight_norm(x61, w_ib61, n);

			Z31::weight(x31, w_ib31, n);
			for (size_t m = n / 4, s = 1; m >= 1; m /= 2, s *= 2) GF31::forward2<1>(x31, m, s);
			GF31::sqr(x31, &w31[n / 4], n);
			for (size_t m = 1, s = n / 4; s >= 1; m *= 2, s /= 2) GF31::backward2<1>(x31, m, s);
			Z31::unweight_norm(x31, w_ib31, n);
		}
		else
		{
			// weighted convolution, radix-2 + radix-3 transforms

			Z61::weight(x61, w_ib61, n);
			for (size_t m = n / 4 / 3, s = 1; m >= 1; m /= 2, s *= 2) GF61::forward2<3>(x61, m, s);
			GF61::sqr3(x61, n);
			for (size_t m = 1, s = n / 4 / 3; s >= 1; m *= 2, s /= 2) GF61::backward2<3>(x61, m, s);
			Z61::unweight_norm(x61, w_ib61, n);

			Z31::weight(x31, w_ib31, n);
			for (size_t m = n / 4 / 3, s = 1; m >= 1; m /= 2, s *= 2) GF31::forward2<3>(x31, m, s);
			GF31::sqr3(x31, n);
			for (size_t m = 1, s = n / 4 / 3; s >= 1; m *= 2, s /= 2) GF31::backward2<3>(x31, m, s);
			Z31::unweight_norm(x31, w_ib31, n);
		}

		// carry propagation
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
			if (x61[k].get() != 0u) return false;
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

// #define	CHECK_MERSENNE_PRIMES	true

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
