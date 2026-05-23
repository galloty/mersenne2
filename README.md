# mersenne2
Search for (large) Mersenne primes using (small) Mersenne primes

## About

Efficient multiplication modulo a Mersenne number is evaluated using an [Irrational Base Discrete Weighted Transform](https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1185244-1/).  
A number theoretic weighted transform can be used. Nick Craig-Wood ([IOCCC 2012 Entry](https://github.com/ncw/ioccc2012/)) implemented it over the field &Zopf;/*p*&Zopf; with *p*&nbsp;=&nbsp;2<sup>64</sup>&nbsp;-&nbsp;2<sup>32</sup>&nbsp;+&nbsp;1.  

**mersenne2** is an implementation of a (large) Mersenne IBDWT over a set of finite fields GF(*M*<sub>*q*</sub><sup>2</sup>), where *M*<sub>*q*</sub> is itself a (small) Mersenne prime. A root of unity of order 3&nbsp;&middot;&nbsp;2<sup>*q*+1</sup> exists in GF(*M*<sub>*q*</sub><sup>2</sup>). Hence, the transform length can be any divisor of 3&nbsp;&middot;&nbsp;2<sup>*q*+1</sup>.

Let *n* be the length of the transform. Over the field &Zopf;/*p*&Zopf;, the NTT is a radix-2 transform, the number of operations is *K*&nbsp;&middot;&nbsp;*n*&nbsp;&middot;&nbsp;log<sub>2</sub>&middot;*n*. Over the field GF(*M*<sub>*q*</sub><sup>2</sup>), the NTT is a radix-4 transform of half length then the number of operations is 3/2&nbsp;&middot;&nbsp;*K*&nbsp;&middot;&nbsp;*n*&nbsp;&middot;&nbsp;log<sub>2</sub>&middot;*n*.  

*M*<sub>61</sub> = 2<sup>61</sup>&nbsp;-&nbsp;1 can be compared to 2<sup>64</sup>&nbsp;-&nbsp;2<sup>32</sup>&nbsp;+&nbsp;1. The number of modular operations is +50% but the modular reduction is easier to calculate. With GF(*M*<sub>*q*</sub><sup>2</sup>), the weights are powers of two: weighting inputs and unweighting outputs are shift operations.  

The important point is that outputs are elements of &Zopf;/*M*<sub>*q*</sub>&Zopf;. A residue number system of two or more Mersenne primes can be used and the solution is obtained by the Chinese Remainder Theorem. With *M*<sub>61</sub> and *M*<sub>31</sub>, the outputs are 92-bit integers. It is more efficient than the single Mersenne prime *M*<sub>89</sub>. GPU are 32-bit processors and an implementation based on (*M*<sub>61</sub>;&nbsp;*M*<sub>31</sub>) is faster than only *M*<sub>61</sub>, if the tested Mersenne number is large enough. The reduction in the size of the transform is a higher gain than the extra 32-bit operations.  
The combination of the two Mersenne primes *M*<sub>61</sub> and *M*<sub>31</sub> is now faster than 2<sup>64</sup>&nbsp;-&nbsp;2<sup>32</sup>&nbsp;+&nbsp;1.  

The Lucas–Lehmer primality test is implemented in C++.  
The Mersenne IBDWT is computed over the two fields GF(*M*<sub>61</sub><sup>2</sup>) and GF(*M*<sub>31</sub><sup>2</sup>).  
The length of the transform is a power of two or three times a power of two.  


## Build

The compiler must support 128-bit literal values. 64-bit version of GCC and Clang support the built-in __uint128_t type.  
The code has been validated using gcc 16 and clang 22.  
