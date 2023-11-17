#!/use/bin/env python3

from argparse import ArgumentParser, ArgumentTypeError

# Sage modules
import sage.all
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.calculus.functional import derivative
from sage.sets.primes import Primes
from sage.arith.misc import divisors
from sage.misc.latex import latex

# Silence warnings about slow implementations in the pre-processing step e.g.
#     verbose 0 (4176: multi_polynomial_ideal.py, groebner_basis)
#     Warning: falling back to very slow toy implementation.
try:
    from sage.misc.verbose import set_verbose
    set_verbose(-1)
except ImportError:
    try:
        from sage.misc.misc import set_verbose
        set_verbose(-1)
    except ImportError:
        pass

(t, ) = PolynomialRing(ZZ, names=('t',))._first_ngens(1)
S = PolynomialRing(ZZ, names=('Q', 'R'))
(Q, R, ) = S._first_ngens(2)


def reduce(poly, n, lift=True):
    # Assumes poly is an element of PolynomialRing(ZZ, names=('t'))
    R = PolynomialRing(ZZ, names=('t'))
    T = R.change_ring(Zmod(n))
    (t, ) = T._first_ngens(1)
    S = T.quotient(t**n)
    if lift:
        return S(poly).lift()
    else:
        return S(poly)

    return poly % n


def compute_weight(p, check=True):
    if p == 0:
        return 0

    if check:
        weights = [4*Q_power + 6*R_power for (Q_power, R_power) in p.exponents()]
        if all([weight == weights[0] for weight in weights]):
            return weights[0]
        print(f'{p} is not homogeneous!')
        exit(1)
    else:
        Q_power, R_power = p.exponents()[0]
        return 4*Q_power + 6*R_power


def init_poly(p, order, OR=False):
    (t, ) = PolynomialRing(ZZ, names=('t'))._first_ngens(1)

    k = compute_weight(p)

    if OR:
        if k != 12:
            print("Error. O'Sullivan-Risager variant only for forms of weight 12")
            print(f"{f} is of weight {k}")
            exit(1)
        return 1

    if p == 0:
        return 0*t

    # # Sage implementation detail
    # # In theory, an approach like this would be ideal
    #
    # R1 = LaurentSeriesRing(QQ, names=('q')); (q, ) = R1._first_ngens(1)
    # R2 = LaurentSeriesRing(R1, names=('r')); (r, ) = R2._first_ngens(1)
    # R3 = PolynomialRing(R3, names=('w')); (w, ) = R3._first_ngens(1)
    #
    # if order == 3:
    #     S = R3.quotient(q**12 * r**(-8) - w)
    #     # ...
    # elif order == 2:
    #     S = R3.quotient(r**12 * q**(-18) - w)
    #     # ...
    #
    # # However, we have the following issue
    # print(S(q*r).lift())
    # q*r    # Here we want w!
    #
    # # So we are forced to use an approach like this (e.g. for order = 2)
    # R = LaurentSeriesRing(QQ, names=('q')); (q, ) = R._first_ngens(1)
    # RR = PolynomialRing(R, names=('r', 'w')); (r, w, ) = RR._first_ngens(2)
    # S = RR.quotient(q*r-w)
    # print(S(q*r).lift())
    # w

    if order == 3:
        # Point = \rho
        # R(\rho) != 0, so we can invert R:
        #     Q^a R^b = R^{k/6} (Q R^{-2/3})^a
        # Sage does not support fractional powers for LaurentSeriesRing
        # Think: r = R^{1/12}
        R = LaurentSeriesRing(QQ, names=('r'))
        (r, ) = R._first_ngens(1)
        # Think: q = Q^{1/12}
        # w is a helper ariable used for the quotient
        RR = PolynomialRing(R, names=('q', 'w'))
        (q, w, ) = RR._first_ngens(2)
        # We want w = QR^{-2/3} = q^{12}r^{-8}
        S = RR.quotient(q**12 * r**(-8) - w)

        s = 0*q + 0*r
        for coef, (Q_power, R_power) in zip(p.coefficients(), p.exponents()):
            s += coef * q**(12*Q_power) * r**(12*R_power)
        s *= r**(-2*k)
        s = S(s).lift()

    elif order == 2:
        # Point = i
        # Q(i) != 0, so we can invert Q:
        #     Q^a R^b = Q^{k/4} (Q^{-3/2} R)^b
        # Sage does not support fractional powers for LaurentSeriesRing
        # Think: q = Q^{1/12}
        R = LaurentSeriesRing(QQ, names=('q'))
        (q, ) = R._first_ngens(1)
        # Think: r = R^{1/12}
        # w is a helper variable used for the quotient
        RR = PolynomialRing(R, names=('r', 'w'))
        (r, w, ) = RR._first_ngens(2)
        S = RR.quotient(r**12*q**(-18) - w)

        s = 0*q + 0*r
        for coef, (Q_power, R_power) in zip(p.coefficients(), p.exponents()):
            s += coef * q**(12*Q_power) * r**(12*R_power)
        s *= q**(-3*k)
        s = S(s).lift()

    p_init = 0*t
    # s is a polynomial in RR (it was sent to S and the lifted back to RR)
    # Therefore, although it is expressed only in the variable w, it is
    # formally in r and w (See the definition of RR)
    # So we throw away the first in every pair using (_, exp)
    for coef, (_, exp) in zip(s.coefficients(), s.exponents()):
        p_init += int(coef.subs(r=0, w=0, q=0))*t**exp

    return p_init


class Recursion:
    def vprint(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)

    def _reduce(self, poly): return reduce(poly, self.d) if self.d else poly

    def __init__(self, f, order, d=None, verbose=True, OR=False):
        self.n = 0
        self.order = order
        self.d = d
        self.verbose = verbose
        self.OR = OR
        self.k = compute_weight(f)

        if OR and self.k != 12:
            print("Error. O'Sullivan-Risager variant only for forms of weight 12")
            print(f"{f} is of weight {self.k}")
            exit(1)

        if self.OR:
            try:
                self.f = S(f/(Q**3 - R**2))
            except TypeError:
                print(f"Error. When using the OR-method, {f} must be given"
                      f"as explicit multiple of the discriminant Q^3-R^2")
                exit(1)
        else:
            self.f = f

        self.p0 = self._reduce(init_poly(self.f, order))
        theta_p0 = -(4*R*derivative(self.f, Q) + 6*Q**2*derivative(self.f, R))
        self.p1 = self._reduce(init_poly(theta_p0, order))

    def __iter__(self): return self

    def __next__(self):
        if self.n == 0:
            self.n += 1
            return self._reduce(self.p0)
        elif self.n == 1:
            self.n += 1
            return self._reduce(self.p1)

        if self.order == 2 and self.OR:
            p2 = self._reduce(6*(t**2-1)*derivative(self.p1)
                              - 2*(self.n-1)*t*self.p1
                              - (self.n-1)*(self.n+10)*self.p0)
        elif self.order == 2 and not self.OR:
            p2 = self._reduce(6*(t**2-1)*derivative(self.p1)
                              - (self.k+2*self.n-2)*t*self.p1
                              - (self.n-1)*(self.n+10)*self.p0)
        elif self.order == 3 and self.OR:
            p2 = self._reduce(4*(t**3-1)*derivative(self.p1)
                              - 2*(self.n-1)*t**2*self.p1
                              - (self.n-1)*(self.n+10)*t*self.p0)
        elif self.order == 3 and not self.OR:
            p2 = self._reduce(4*(t**3-1)*derivative(self.p1)
                              - (self.k+2*self.n-2)*t**2*self.p1
                              - (self.n-1)*(self.n+self.k-2)*t*self.p0)

        self.n += 1
        self.p0, self.p1 = self.p1, p2
        return self._reduce(self.p1)


def minimal_period(seq, period):
    # Given a periodic (from the start) sequence and a period, find the
    # minimal period
    for min_period in divisors(period):
        # Use a generator expression for lazy eval
        if all((seq[:min_period] == seq[min_period*k:min_period*(k+1)]
                for k in range(1, period//min_period)
                )):
            break
    return min_period


def period_offset(R, r, d):
    # Compute the first r terms
    S = [next(R) for _ in range(r)]
    P = 0
    while not P:
        # Append the next value in the sequence
        S += [next(R)]
        # Move a window of r elements over the sequence to verify condition (4.2)
        for B in range((len(S)-r) % d, len(S)-r, d):
            if S[B:B+r] == S[len(S)-r:len(S)]:
                # Condition (4.2) is met with M = B, N = len(S) - r
                # Sequence is `N - M` periodic with offset `B`
                P = len(S) - r - B
                break
    # Compute the minimal period `p`
    for p in divisors(P):
        if all([S[B:B+p] == S[B+p*k:B+p*(k+1)] for k in range(1, P//p)]):
            break
    # Compute another `p` terms to ensure `S` contains at least two periods
    S += [next(R) for _ in range(p)]
    # Compute minimal offset
    for b in range(B+1):
        if S[b:b+p] == S[b+p:b+2*p]:
            break
    return p, b, S


def pprint(polys, d, period, offset, rd, LaTeX=False):
    # Pretty print the polynomials

    poly_strs = [f"p_{{{n}}}(t) = {reduce(poly, d)}"
                 for n, poly in enumerate(polys[:offset+period+rd])]
    max_len = max([len(poly_str) for poly_str in poly_strs])

    if LaTeX:
        print("\\noindent")
        print(f"Modulo {d}, the polynomials have period {period} and offset {offset}")
        print()

        if offset:
            print("\\noindent")
            print(f"The first {offset} polynomials are")
            print(f"% The first {offset} polynomials are " + "{{{{{")
            print()
            print('\\begin{multicols}{2}')
            for n, poly in enumerate(polys[:offset]):
                print(f"  \\noindent $p_{{{n}}} \\equiv {latex(poly)}$\n")
            print('\\end{multicols}')
            print("% }}}}}")

        print()
        print("\\noindent")
        print(f"The {period} repeating polynomials are")
        print(f"% The {period} repeating polynomials are " + "{{{{{")

        print('\\begin{multicols}{2}')
        for n, poly in enumerate(polys[offset:offset+period]):
            print(f"  \\noindent $p_{{{n+offset}}} \\equiv {latex(poly)}$\n")
        print('\\end{multicols}')
        print("% }}}}}")
        print()

    else:
        if offset:
            print()
            print(f"The first {offset} non-periodic polynomials")
            print()
            for poly in poly_strs[:offset]:
                print(f"  {poly.ljust(max_len)} (mod {d})")
            print()

        print(f"Here is a repeating period of {period} polynomials")
        print("This sequence repeats forever")
        print()

        for poly in poly_strs[offset:offset+period]:
            print(f"  {poly.ljust(max_len)} (mod {d})")

        for n, poly in enumerate(poly_strs[offset+period:offset+period+rd]):
            print(f"  {poly.ljust(max_len)} (mod {d}) = {poly_strs[n+offset]}")
        print("  ...")


def pprint_p0(polys, d, order, period, offset, LaTeX=False):
    # Pretty print the polynomials evaluated at 0

    if not LaTeX:
        print("Now finding periodic behaviour of p_n(t=0)")

    # Since the polys are `period` periodic (after offset `offset`), the same
    # must be true for evaluating at 0, however p(0) may have a shorter period
    values = [int(poly.subs(t=0)) for poly in polys[offset:offset+period]]
    s_period = minimal_period(values, period)

    if LaTeX:
        print("\\noindent")
        print(f"Values of $p_{{n}}(0) \\pmod{{{d}}}$")
        print(f"% Values of $p_{{n}}(0) \\pmod{{{d}}}$" + "{{{{{")
        print('\\begin{multicols}{3}')
        for n, poly in enumerate(polys[offset:offset+s_period]):
            if poly.subs(t=0) != 0:
                print(f"  \\noindent $p_{{{s_period}n + {(n+offset)%s_period}}}(0)"
                      f" \\equiv {poly.subs(t=0)} \\pmod{{{d}}}$")
                print()
        print('\\end{multicols}')
        print("% }}}}}")
        print()

    else:
        print(f"The values p_{{n}}(0) are {s_period} periodic")
        print(f"That is, the non-trivial values p_{{n}}(0) are {s_period//order} periodic:")

        for n, poly in enumerate(polys[offset:offset+s_period]):
            if poly.subs(t=0) != 0:
                print(f"  p_{{{s_period}*n + {(n+offset)%s_period}}}(0)"
                      f" = {poly.subs(t=0)} (mod {d})")


def compute(f, order, candidates=[5], max_candidates=1, verbose=True, OR=False, LaTeX=False):
    def vprint(*args, **kwargs):
        if verbose and not LaTeX:
            print(*args, **kwargs)

    k = compute_weight(f)

    return_str = []
    successful_candidates = 0

    for d in candidates:
        # Initialise recursion
        R = Recursion(f, order, d=d, verbose=verbose, OR=OR)
        # Get periodicity, offset, and the computed polynomials
        # `comp` includes (at least) two full periods
        period, offset, comp = period_offset(R, 2, d)
        if comp[offset:offset+period] != comp[offset+period:offset+2*period]:
            print("Failed period computation")

        msg = ""
        msg += f"weight = {k}"
        msg += f", fn = {f}"
        msg += f", order = {order}"
        msg += f", candidate = {d:>4}"
        msg += f", offset = {offset:>4}"
        msg += f", period = {period:>4}"
        msg += f" (took {len(comp)} elements)"

        zero = False
        for n, poly in enumerate(comp[offset + abs(k//2 - offset) % order::order]):
            if poly.subs(t=0) == 0:
                nth = n*order + offset + abs(k//2 - offset) % order
                msg += ", Not all non-trivial Fourier coefficients are nonzero"
                msg += f", e.g. the {nth}th poly p_{{{nth}}}(t) = {poly} (mod {d})"
                zero = True
                break

        if not zero:
            # vprint(f"All non-trivial Fourier coefficients are non-zero modulo {d}")
            successful_candidates += 1
            msg += f", All non-trivial Fourier coefficients are non-zero modulo {d}"

        vprint(msg)
        return_str += [msg]

        if verbose and not zero:
            pprint(comp, d, period, offset, 2, LaTeX=LaTeX)
            pprint_p0(comp, d, order, period, offset, LaTeX=LaTeX)

        if successful_candidates == max_candidates:
            break

    return return_str


def qr(polynomial):
    """Verifies that the polynomial given is one of Q, R"""
    S = PolynomialRing(ZZ, names=('Q', 'R'))
    (Q, R, ) = S._first_ngens(2)
    try:
        S(polynomial)
    except TypeError:
        raise ArgumentTypeError(f"'{polynomial}' is not a polynomial in Q, R")
    return S(polynomial)


def main():
    parser = ArgumentParser(description='Prints polynomials for computing the Taylor Expansion'
                            ' of a modular form at a given point.'
                            ' The periodicity of the polynomials is computed'
                            )
    parser.add_argument('--OR', action='store_true',
                        help="Use O'Sullivan-Risager's method."
                             " Only works when function is given as explicit"
                             " multiple of Q^3 - R^2."
                        )
    parser.add_argument('--latex', action='store_true',
                        help="Pretty print output in LaTeX compatible format"
                        )
    parser.add_argument('f', type=qr,
                        help="Function to expand."
                             " Expressed as a polynomial in R, Q"
                        )
    parser.add_argument('order', type=int, choices=[2, 3],
                        help="Order of the point to expand at."
                             " 2: Expands at i. 3: Expands at \\rho"
                        )
    args = parser.parse_args()

    # Exclude 2 and 3 for now (2, 3 are the only factors of 12)
    P = Primes().unrank_range(2, 1000)
    (t, ) = PolynomialRing(ZZ, names=('t',))._first_ngens(1)
    S = PolynomialRing(ZZ, names=('Q', 'R'))
    (Q, R, ) = S._first_ngens(2)

    compute(args.f,
            args.order,
            candidates=P,
            OR=args.OR,
            verbose=True,
            LaTeX=args.latex)


if __name__ == "__main__":
    main()
