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
