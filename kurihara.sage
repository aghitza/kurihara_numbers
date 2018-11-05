LOGS = dict()


def kurihara_number(E, p, n):
    """
    Return the Kurihara number for the elliptic curve E with respect to
    the prime p and the product of Kolyvagin primes n.

    Note that the Kurihara number itself is not well-defined (it depends
    on choices of multiplicative generators), but its vanishing and
    non-vanishing are well-defined.
    """
    f = E.modular_symbol()
    ell1, ell2 = n.prime_divisors()
    K1 = GF(ell1)
    g1 = K1.multiplicative_generator()
    if not ell1 in LOGS:
        LOGS[ell1] = precompute_logs(g1)
    K2 = GF(ell2)
    g2 = K2.multiplicative_generator()
    if not ell2 in LOGS:
        LOGS[ell2] = precompute_logs(g2)
    K = GF(p)
    S = K(0)
    for a in range(n):
        if gcd(a, n) == 1:
            aa1 = K1(a)
            aa2 = K2(a)
            S = S + K(LOGS[ell1][aa1]) * K(LOGS[ell2][aa2]) * K(f(a/n))
    return S


def precompute_logs(g):
    K = g.parent()
    logs = dict()
    for a in K:
        if a != 0:
            logs[a] = a.log(g)
    return logs


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
