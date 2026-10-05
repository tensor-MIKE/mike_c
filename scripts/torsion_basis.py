from sage.all import ZZ, GF, EllipticCurve, parallel

def even_torsion_basis_E0(E0, f):
    """
    For the case when A = 0 we can't use the entangled basis algorithm
    so we do something "stupid" to simply get something canonical
    """
    assert E0.a_invariants() == (0, 0, 0, 1, 0)

    Fp2 = E0.base_ring()
    p = Fp2.characteristic()

    def points_order_two_f():
        """
        Compute a point P of order 2^f with x(P) = 1 + i*x_im
        """
        x_im = 0
        while True:
            x_im += 1
            x = Fp2([1, x_im])
            if not E0.is_x_coord(x):
                continue
            # compares a+bi <= c+di iff (a,b) <= (c,d) as tuples, where integers
            # modulo p are compared via their minimal non-negative representatives
            P = min(E0.lift_x(x, all=True), key = lambda pt: list(pt.y()))
            P.set_order(multiple=p+1)
            if P.order() % (1 << f) == 0:
                P *= P.order() // (1 << f)
                P.set_order(1 << f)
                yield P

    pts = points_order_two_f()
    P = next(pts)
    for Q in pts:
        # Q is picked to be in E[2^f] AND we must ensure that
        # <P, Q> form a basis, which is the same as e(P, Q) having
        # full order 1 << f.
        e = P.weil_pairing(Q, 1 << f)
        if e ** (1 << f - 1) == -1:
            break

    # Finally we want to make sure Q is above (0, 0)
    P2 = (1 << f - 1) * P
    Q2 = (1 << f - 1) * Q
    if Q2 == E0(0, 0):
        pass
    elif P2 == E0(0, 0):
        P, Q = Q, P
    else:
        Q += P

    assert P.weil_pairing(Q, 1 << f) ** (1 << f - 1) == -1
    assert (1 << f - 1) * Q == E0(0, 0)

    return P, Q


def precompute_E0_basis(E0,f):
    """
    Compute the basis E0[2^(f - 1)]
    We use the method of finding a point P with Fp rational coordinates
    and then compute Q = (-x, i*y)
    """
    Fp2 = E0.base_ring()
    assert(E0 == EllipticCurve(Fp2, [-1, 0]))
    p = Fp2.characteristic()
    odd_cofactor = (p + 1) // 2**f
    Fp = GF(p)
    assert(Fp!=Fp2)
    x = Fp.zero()
    while True:
        # x must be an element of Fp such that (1 - x) is a NQR
        x += Fp.one()
        if (1 - x).is_square():
            continue
        # P must be a point with Fp rational y in E[2^f-1]
        P = odd_cofactor * E0.lift_x(x)
        if P.y() not in Fp:
            continue
        T = 2 ** (f - 2) * P
        if not T:
            continue
        break

    Q = E0(-P.x(), Fp2.gen() * P.y())

    assert not 2 ** (f - 1) * P and not 2 ** (f - 1) * Q
    assert 2 ** (f - 2) * P and 2 ** (f - 2) * Q
    assert not (2 ** (f - 2) * P == 2 ** (f - 2) * Q)

    return P, Q

def EM_basis_E0(E0, f):
    """
    We have the curve and basis:
    E0 : y^2 = x(x + 1)(x - 1), <P, Q> = E[2^(f - 1)]
    with P above (1, 0) and Q above (-1, 0) which means
    the kernel P + [x]Q will be above (0, 0)
    The idea of this map is to compute E0_M in the Montgomery
    model, chosen such that we map Q to be above (0, 0) to ensure
    that regardless of x, P + [x]Q is never above (0, 0).
    To do this, we pick alpha, a root of  x(x + 1)(x - 1) to be 0, 1, or -1.
    We then compute s = 1 / (3*alpha^2 - 1) such that
    EM : B * y^2 = x^3 + Ax^2 + x; A = 3 alpha s, B = s
    (x, y) -> (s(x - alpha), sy)
    We see that to ensure that Q is mapped above (0, 0) we must select
    alpha = 1.
    We also want B = 1, so we need to have another factor of u = s.sqrt()
    when we compute y.
    """
    Fp2 = E0.base_ring()
    
    alpha = 1  # Pick this root to send Q to be above zero
    s = ~Fp2(3 * alpha**2 - 1).sqrt()
    u = s.sqrt()
    A = -3 * alpha * s
    EM = EllipticCurve(Fp2, [0, A, 0, 1, 0])
    def iso(P):
        x, y = P.xy()
        x_new = (x + alpha) * s
        y_new = y * (s * u)  # Ensure B = 1
        return EM(x_new, y_new)
    # Compute the new basis from the special E0 basis
    EM_basis = [iso(P).division_points(2)[0] for P in precompute_E0_basis(E0,f)]
    # Ensure the mapping has worked as expected
    assert 2 ** (f - 1) * EM_basis[0] != EM(0, 0)
    assert 2 ** (f - 1) * EM_basis[1] == EM(0, 0)
    return EM, EM_basis[0], EM_basis[1]
