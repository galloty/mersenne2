# mersenne2
Search for (large) Mersenne primes using (small) Mersenne primes

## About

Efficient multiplication modulo a Mersenne number is evaluated using an [Irrational Base Discrete Weighted Transform](https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1185244-1/).  
A number theoretic weighted transform can be used. Nick Craig-Wood ([IOCCC 2012 Entry](https://github.com/ncw/ioccc2012/)) implemented it over the field Z/*p*Z with *p* = 2<sup>64</sup>&nbsp;-&nbsp;2<sup>32</sup>&nbsp;+&nbsp;1.  

**mersenne2** is an implementation of (large) Mersenne IBDWT over the finite field GF(*M*<sub>*q*</sub><sup>2</sup>), where *M*<sub>*q*</sub> is itself a (small) Mersenne prime. We have 2<sup>q</sup> = 1 (mod&nbsp;*M*<sub>*q*</sub>) then a *n*<sup>th</sup> root of two always exists if *n* is a power of two. A root of unity of order 2<sup>*q*+1</sup> exists in GF(*M*<sub>*q*</sub><sup>2</sup>). Hence, the transform length can be any divisor of 2<sup>*q*+1</sup>.

Let *n* be the length of the transform. Over the field Z/*p*Z, the NTT is a radix-2 transform, the number of operations is *K*&nbsp;&middot;&nbsp;*n*&nbsp;&middot;&nbsp;log<sub>2</sub>&middot;*n*. Over the field GF(*M*<sub>*q*</sub><sup>2</sup>), the NTT is a radix-4 transform of half length then the number of operations is 3/2&nbsp;&middot;&nbsp;*K*&nbsp;&middot;&nbsp;*n*&nbsp;&middot;&nbsp;log<sub>2</sub>&middot;*n*.  

*M*<sub>61</sub> = 2<sup>61</sup>&nbsp;-&nbsp;1 can be compared to 2<sup>64</sup>&nbsp;-&nbsp;2<sup>32</sup>&nbsp;+&nbsp;1. The number of modular operations is +50% but the modular reduction is easier to calculate. With GF(*M*<sub>*q*</sub><sup>2</sup>), the weights are powers of two: weighting inputs and unweighting outputs are shift operations.  

An important point: outputs are elements of Z/*q*Z. Two or more Mersenne primes can be used. We have a residue number system and the solution is obtained by the Chinese Remainder Theorem. With *M*<sub>61</sub> and *M*<sub>31</sub>, the outputs are 92-bit integers. It is more efficient then the single Mersenne prime *M*<sub>89</sub>. GPU are 32-bit processors and an implementation based on (*M*<sub>61</sub>;&nbsp;*M*<sub>31</sub>) is faster than only *M*<sub>61</sub> if the tested Mersenne number is large enough. The reduction in the size of the transform is a higher gain than the extra 32-bit operations.
