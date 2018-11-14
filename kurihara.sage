# cache Kurihara numbers
KN = dict()


def kurihara_graph(E, p, plst, verbose=False):
    res = []
    from sage.combinat.subset import SubsetsSorted
    for s in SubsetsSorted(plst):
        n = prod(s)
        if verbose:
            print("Computing Kurihara number for n=%s (%s)" % (n, s))
        k = kurihara_number(E, p, n, verbose=verbose)
        if k == 0:
            res.append(0)
        else:
            res.append(1)
    return tuple(res)


def tikz_from_kurihara_graph(E, p, plst, kg, verbose=False):
    pass


def kurihara_number(E, p, n, verbose=False):
    """
    Return the Kurihara number for the elliptic curve E with respect to
    the prime p and the product of Kolyvagin primes n.

    Note that the Kurihara number itself is not well-defined (it depends
    on choices of multiplicative generators), but its vanishing and
    non-vanishing are well-defined.
    """
    if n in KN:
        return KN[n]
    LOGS = dict()
    from sage.libs.eclib.newforms import ECModularSymbol
    f = ECModularSymbol(E, sign=int(1))
    # we're using the eclib implementation of modular symbols directly,
    # and the normalisation is different than Sage's if the curve E has
    # negative discriminant, hence the scaling factor
    scale = 1
    if E.discriminant() < 0:
        scale = 2
    ells = n.prime_divisors()
    Kell = dict()
    if verbose:
        print("precomputing discrete logs...")
    for ell in ells:
        Kell[ell] = GF(ell)
        if not ell in LOGS:
            LOGS[ell] = precompute_logs(ell)
    if verbose:
        print("...done")
    K = GF(p)
    S = K(0)
    for a in xsrange(1, n):
        if gcd(a, n) == 1:
            if verbose:
                print(RDF(a/n)*100)
            mult = K(f(a/n))
            for ell in ells:
                mult *= K(LOGS[ell][Kell[ell](a)])
            S += mult
    res = S/scale
    KN[n] = res
    return res


def precompute_logs(ell):
    K = GF(ell)
    g = K.multiplicative_generator()
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
