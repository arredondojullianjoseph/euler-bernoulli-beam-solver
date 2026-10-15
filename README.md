# Euler-Bernoulli Beam Solver

One-element cubic Hermite cantilever. Galerkin $k$ vs Hibbeler $PL^3/(3EI)$. Measured relative error $2.6\times 10^{-16}$.

**Next Up:** simply-supported load, then a mesh check.

**Author:** Jullian J. Arredondo (jjarredondo@liberty.edu)

## Model

$$
EI w'''' = q,\qquad k_{ij} = \int EI\, N_i'' N_j''\,\mathrm{d}x,\qquad w_{\mathrm{tip}} = \frac{PL^3}{3EI}.
$$

Wall at $x=0$. Tip load $P$. Paper numbers: $L=1\,\mathrm{m}$, $E=200\,\mathrm{GPa}$, $I=10^{-6}\,\mathrm{m}^4$, $P=1000\,\mathrm{N}$.

## Files

- `src/closed_form.py` — $L,E,I,P$ and Hibbeler
- `src/stiffness.py` — 4×4 $k$
- `src/beam.py` — wall, load, solve

## Run

```bash
pip install -r requirements.txt
python src/beam.py
```

| w_tip (1 element) | Hibbeler PL^3/(3EI) | relative error |
| --- | --- | --- |
| 1.666666666667e-03 m | 1.666666666667e-03 m | 2.6e-16 |

Not Timoshenko, not a frame, not commercial FEA.

Hibbeler RC. *Mechanics of Materials*.
