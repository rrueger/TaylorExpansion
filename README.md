TaylorExpansion
--------------

A python (sage) module that computes the Taylor coefficients (up to scaling) of
a modular form at the CM points $i$ and $\rho$.

It also supplies the command-line program `taylor-expansion` with the following
usage.

```
usage: taylor-expansion [-h] [--OR] [--latex] f {2,3}

Prints polynomials for computing the Taylor Expansion of a modular form at a
given point. The periodicity of the polynomials is computed

positional arguments:
  f           Function to expand. Expressed as a polynomial in R, Q
  {2,3}       Order of the point to expand at. 2: Expands at i. 3: Expands at
              \rho

options:
  -h, --help  show this help message and exit
  --OR        Use O'Sullivan-Risager's method. Only works when function is
              given as explicit multiple of Q^3 - R^2.
  --latex     Pretty print output in LaTeX compatible format
```

Install for usage with `python` with

    pip install https://github.com/rrueger/TaylorExpansion/raw/main/dist/sage_taylorexpansion-0.9.0-py3-none-any.whl

Sage uses its own package hierarchy to match the version of `python` it is
shipped with. To install for usage with `sage` call

    sage -pip install https://github.com/rrueger/TaylorExpansion/raw/main/dist/sage_taylorexpansion-0.9.0-py3-none-any.whl

*Note* Sage is not officially available as a python module on PyPi, therefore it
is not listed as a dependency in this module and will **not** be automatically
installed alongside `TaylorExpansion`. Sage must be installed on the system
independently, before `TaylorExpansion` can be installed with `sage -pip` as
above.

### Example

    taylor-expansion 'Q^3 - R^2' 2


```
weight = 12, fn = Q^3 - R^2, order = 2, candidate =    5, offset =    0, period =   20 (took 42 elements), All non-trivial Fourier coefficients are non-zero modulo 5
Here is a repeating period of 20 polynomials
This sequence repeats forever

  p_{0}(t) = 4*t^2 + 1          (mod 5)
  p_{1}(t) = 0                  (mod 5)
  p_{2}(t) = 2*t^2 + 3          (mod 5)
  p_{3}(t) = 2*t^3 + 3*t        (mod 5)
  p_{4}(t) = 4*t^2 + 1          (mod 5)
  p_{5}(t) = 3*t^3 + 2*t        (mod 5)
  p_{6}(t) = 3*t^4 + 4*t^2 + 3  (mod 5)
  p_{7}(t) = 4*t^3 + t          (mod 5)
  p_{8}(t) = 4*t^2 + 1          (mod 5)
  p_{9}(t) = 3*t^3 + 2*t        (mod 5)
  p_{10}(t) = 4*t^4 + 3*t^2 + 3 (mod 5)
  p_{11}(t) = 4*t^3 + 3*t       (mod 5)
  p_{12}(t) = 3*t^4 + 3*t^2 + 1 (mod 5)
  p_{13}(t) = 2*t^3             (mod 5)
  p_{14}(t) = 4*t^4 + 3*t^2 + 3 (mod 5)
  p_{15}(t) = 4*t               (mod 5)
  p_{16}(t) = t^2 + 1           (mod 5)
  p_{17}(t) = 3*t^3 + t         (mod 5)
  p_{18}(t) = t^4 + 3           (mod 5)
  p_{19}(t) = 4*t               (mod 5)
  p_{20}(t) = 4*t^2 + 1         (mod 5) = p_{0}(t) = 4*t^2 + 1
  p_{21}(t) = 0                 (mod 5) = p_{1}(t) = 0
  ...
Now finding periodic behaviour of p_n(t=0)
The values p_{n}(0) are 4 periodic
That is, the non-trivial values p_{n}(0) are 2 periodic:
  p_{4*n + 0}(0) = 1 (mod 5)
  p_{4*n + 2}(0) = 3 (mod 5)
```

This can be verified to be the same as the result Proposition 5.1 from the Paper
_Non-vanishing of Taylor Coefficients and Poincare series_ by O'Sullivan and
Risager.
