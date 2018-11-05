def kurihara_number(E, p, n):
    """
    Return the Kurihara number for the elliptic curve E with respect to
    the prime p and the product of Kolyvagin primes n.

    Note that the Kurihara number itself is not well-defined (it depends
    on choices of multiplicative generators), but its vanishing and
    non-vanishing are well-defined.
    """
    f = E.modular_symbol()
    R1 = Integers(factor(n)[0][0])
    g1 = R1.multiplicative_generator()
    R2 = Integers(factor(n)[1][0])
    g2 = R2.multiplicative_generator()
    S = mod(0, p)
    for a in range(n):
        if gcd(a, n) == 1:
            aa1 = mod(a, factor(n)[0][0])
            aa2 = mod(a, factor(n)[1][0])
            S = S + mod(aa1.log(g1), p) * mod(aa2.log(g2), p) * mod(f(a/n), p)
    return S


def kolyvagin_primes(E, p, bound):
    """
    Return the list of Kolyvagin primes less than bound for the elliptic
    curve E with respect to the prime p.
    """
    lst = []
    for ell in prime_range(bound):
        if mod(ell-1, p) == 0:
            if mod(E.ap(ell) - ell - 1, p) == 0:
                lst.append(ell)
    return lst


def check_prime(E, p):
    """
    Return True if p satisfies Assumption 1.1.
    """
    if p == 2:
        return False
    N = E.conductor()
    if N % p == 0:
        return False
    if E.tamagawa_product() % p == 0:
        return False
    for ell in N.prime_divisors():
        if E.has_nonsplit_multiplicative_reduction(ell):
            if (ell+1) % p == 0:
                return False
        if E.has_split_multiplicative_reduction(ell):
            if (ell-1) % p == 0:
                return False
    return True
