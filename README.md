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
- Simply-supported: pins on $w$ at both ends, uniform $q$. Hibbeler $w_{\mathrm{mid}}=5qL^4/(384EI)$. Exact field $w(x)=(qx/24EI)(L^3-2Lx^2+x^3)$.

## Implementation

Same cubic Hermite $k$ as before. Now $n$ equal pieces of length $h=L/n$. Each node stores $w$ and $\theta$. Uniform $q$ still uses $qh[1/2,\ h/12,\ 1/2,\ -h/12]$ on every piece. Wall: drop $w$ and $\theta$ at $x=0$. Pins: drop only $w$ at both ends.

Nodes are exact for these loads. Error is between nodes. One element at midspan is $qL^4/96EI=4/5$ of Hibbeler. Two elements put a node at $L/2$, so that check is roundoff.

Mesh $2,4,8,16$. $L/4$ is not a node on 2 elements (~1.8% low). On 4/8/16 it is a node, so that error is junk for order. The order number is $\max|w_{\mathrm{FEM}}-w_{\mathrm{exact}}|$ on 101 samples.

## Verification

```bash
pip install -r requirements.txt
python src/beam.py
pytest -s -v
```

| Case | FEM | Hibbeler | relative error |
| --- | --- | --- | --- |
| cantilever $w_{\mathrm{tip}}$ | 1.666666666667e-03 m | 1.666666666667e-03 m | 2.6e-16 |
| SS $w(L/2)$ (2 elements) | 6.510416666667e-05 m | 6.510416666667e-05 m | 2.1e-16 |

$L/4$ and max sample error. Factor is $\mathrm{error}(n)/\mathrm{error}(2n)$ on $\max|dw|$.

| $n$ | $h$ (m) | $w(L/4)$ FEM | $w(L/4)$ exact | rel at $L/4$ | $\max\|dw\|$ (m) | factor |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 0.50 | 4.557291666667e-05 | 4.638671875000e-05 | 1.754e-02 | 8.138e-07 | — |
| 4 | 0.25 | 4.638671875000e-05 | 4.638671875000e-05 | 2.6e-15 | 5.070e-08 | 16.05 |
| 8 | 0.125 | 4.638671875000e-05 | 4.638671875000e-05 | 4.1e-14 | 3.169e-09 | 16.00 |
| 16 | 0.0625 | 4.638671875000e-05 | 4.638671875000e-05 | 4.4e-14 | 1.980e-10 | 16.00 |

Two-element $L/4$ is $1.75\%$ low (about $1.8\%$). After that $L/4$ is a node and the relative error is roundoff. The sample max drops by about $16$ when $h$ is halved.

Plot: `error_vs_h.png` (log-log $\max|dw|$ vs $h$).

## Limitations

- 1D beam only.
- Euler-Bernoulli: no shear; not Timoshenko.
- Linear elastic.
- Static only; no vibration, no time integration.
- Not a frame.
- Not commercial FEA.

## Files

- `src/closed_form.py` — $L,E,I,P,q$, Hibbeler, and the SS field
- `src/stiffness.py` — 4×4 $k$
- `src/beam.py` — cantilever, $n$-element simply-supported, $w(x)$ from $N$, mesh sweep
- `tests/test_beam.py` — tip, midspan, and $L/4$ / mesh-order checks
- `error_vs_h.png` — observed-order plot
