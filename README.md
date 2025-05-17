# mersenne2
Search for (large) Mersenne primes using (small) Mersenne primes

## About

Efficient multiplication modulo a Mersenne number is evaluated using an [Irrational Base Discrete Weighted Transform](https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1185244-1/).  
A number theoretic weighted transform can be used. Nick Craig-Wood ([IOCCC 2012 Entry](https://github.com/ncw/ioccc2012/)) implemented it over the field Z/*p*Z with *p* = 2<sup>64</sup> - 2<sup>32</sup> + 1.  
**mersenne2** is an implementation of (large) Mersenne IBDWT over the finite field GF(*M*<sub>*p*</sub><sup>2</sup>), where *M*<sub>*p*</sub> is itself a (small) Mersenne prime.  