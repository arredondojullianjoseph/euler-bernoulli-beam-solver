# Euler-Bernoulli Beam Solver

[![tests](https://github.com/arredondojullianjoseph/euler-bernoulli-beam-solver/actions/workflows/tests.yml/badge.svg)](https://github.com/arredondojullianjoseph/euler-bernoulli-beam-solver/actions/workflows/tests.yml)

Cubic Hermite beam. Galerkin $k$ vs Hibbeler. Cantilever: one element, tip relative error $2.6\times 10^{-16}$. Simply-supported: two equal elements so $L/2$ is a node.

**Author:** Jullian J. Arredondo (jjarredondo@liberty.edu)

## Model

$$
EI w'''' = q,\qquad k_{ij} = \int EI\, N_i'' N_j''\,\mathrm{d}x.
$$

Paper numbers: $L=1\,\mathrm{m}$, $E=200\,\mathrm{GPa}$, $I=10^{-6}\,\mathrm{m}^4$, $P=1000\,\mathrm{N}$, $q=1000\,\mathrm{N/m}$.

- Cantilever: wall at $x=0$, tip $P$. Hibbeler $w_{\mathrm{tip}}=PL^3/(3EI)$.
- Simply-supported: pins on $w$ at both ends, uniform $q$. Two elements. Hibbeler $w_{\mathrm{mid}}=5qL^4/(384EI)$.

## Files

- `src/closed_form.py` — $L,E,I,P,q$ and both Hibbeler formulas
- `src/stiffness.py` — 4×4 $k$
- `src/beam.py` — cantilever, simply-supported, $w(x)$ from $N$
- `tests/test_beam.py` — both closed-form checks

## Run

```bash
pip install -r requirements.txt
python src/beam.py
pytest -s -v
```

| Case | FEM | Hibbeler | relative error |
| --- | --- | --- | --- |
| cantilever $w_{\mathrm{tip}}$ | 1.666666666667e-03 m | 1.666666666667e-03 m | 2.6e-16 |
| SS $w(L/2)$ (2 elements) | 6.510416666667e-05 m | 6.510416666667e-05 m | 2.1e-16 |

One cubic element is exact at the nodes for these loads. Midspan is a node only after two elements; that is why the SS check uses two.
