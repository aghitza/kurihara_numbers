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


def kurihara_number(E, p, n, recalculate=False, verbose=False):
    """
    Return the Kurihara number for the elliptic curve E with respect to
    the prime p and the product of Kolyvagin primes n.

    Note that the Kurihara number itself is not well-defined (it depends
    on choices of multiplicative generators), but its vanishing and
    non-vanishing are well-defined.
    """
    label = E.label() + '_' + str(p) + '__' + '_'.join([str(ell) for ell in n.prime_divisors()])
    fname = 'data/' + label + '.txt'
    if not recalculate and os.path.isfile(fname):
        with open(fname, 'r') as f:
            st = f.readline()
            return GF(p)(st.strip())
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
    L = GF(p)
    if verbose:
        print("precomputing discrete logs...")
    for ell in ells:
        Kell[ell] = GF(ell)
        if not ell in LOGS:
            LOGS[ell] = precompute_logs(L, ell)
    if verbose:
        print("...done")
    S = L(0)
    for a in xsrange(1, n):
        if gcd(a, n) == 1:
            if verbose:
                print(RDF(a/n)*100)
            mult = L(f(a/n))
            for ell in ells:
                mult *= LOGS[ell][Kell[ell](a)]
            S += mult
    res = S/scale
    with open(fname, 'w') as f:
        f.write(str(res) + '\n')
    return res


def precompute_logs(L, ell):
    K = GF(ell)
    g = K.multiplicative_generator()
    logs = dict()
    for a in K:
        if a != 0:
            logs[a] = L(a.log(g))
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
    Return True if the pair (E, p) satisfies Assumptions 1.1 and 1.2.
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
    if E.manin_constant() != 1:
        return False
    rho = E.galois_representation()
    if not rho.is_surjective(p):
        return False
    return True
