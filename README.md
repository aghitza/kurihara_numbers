This repository contains the Sage code for performing the computations described in the Appendix of the paper *Indivisibility of Kato's Euler systems and Kurihara numbers* by Chan-Ho Kim.

The code computes the Kurihara number corresponding to an elliptic curve E, a prime p, and a squarefree product of Kolyvagin primes for the pair (E, p).

The example described in the paper's Appendix can be obtained using Sage as follows (tested with versions 8.6 and 8.7 of Sage):

```
sage: %attach kurihara.sage
sage: E = EllipticCurve('128a1')
sage: check_prime(E, 3)
True
sage: lst = kolyvagin_primes(E, 3, 150); lst
[7, 37, 67, 73, 103]
sage: g = kurihara_graph(E, 3, lst)
sage: g
(0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1,
 0, 0, 1, 0, 0, 0, 0, 0, 1)
```

The code saves the Kurihara numbers that it computes under a subdirectory data/ and reuses them if they are required in further computations.
