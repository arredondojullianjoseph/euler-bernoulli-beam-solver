"""
beam.py

Cantilever: one element vs Hibbeler PL^3 / 3EI.
Simply-supported: n equal elements, uniform q.
"""

from pathlib import Path

import numpy as np

from closed_form import E, I, L, P, Q, cantilever_tip, simply_supported_field, simply_supported_midspan
from stiffness import element_stiffness


def hermite_shape(xi, le):
    """Cubic Hermite N at xi = x/le. Order: w_left, theta_left, w_right, theta_right."""
    return np.array(
        [
            1.0 - 3.0 * xi**2 + 2.0 * xi**3,
            le * (xi - 2.0 * xi**2 + xi**3),
            3.0 * xi**2 - 2.0 * xi**3,
            le * (-(xi**2) + xi**3),
        ]
    )


def evaluate_w(u, x, length):
    """Deflection at x. Find the element, then mix the four local N."""
    n_el = len(u) // 2 - 1
    le = length / n_el
    e = int(x / le)
    if e >= n_el:
        e = n_el - 1
    xi = (x - e * le) / le
    return float(hermite_shape(xi, le) @ u[2 * e : 2 * e + 4])


def cantilever():
    """Wall DOFs 0,1 dropped. P on tip w (slot 2)."""
    k = element_stiffness(E * I, L)
    f = np.zeros(4)
    f[2] = P

    k_free = k[2:, 2:]
    f_free = f[2:]
    u_free = np.linalg.solve(k_free, f_free)

    u = np.zeros(4)
    u[2:] = u_free
    return u


def simply_supported(n_el=2):
    """Pins on w at both ends. Uniform q. n equal elements. Default 2 so L/2 is a node."""
    le = L / n_el
    n_dof = 2 * (n_el + 1)
    K = np.zeros((n_dof, n_dof))
    F = np.zeros(n_dof)
    ke = element_stiffness(E * I, le)
    fe = Q * le * np.array([0.5, le / 12.0, 0.5, -le / 12.0])
    for e in range(n_el):
        slots = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        for a in range(4):
            F[slots[a]] += fe[a]
            for b in range(4):
                K[slots[a], slots[b]] += ke[a, b]
    free = [i for i in range(n_dof) if i not in (0, n_dof - 2)]
    u = np.zeros(n_dof)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    return u


def mesh_sweep():
    """2, 4, 8, 16. Each row: n_el, h, w_fem(L/4), w_exact(L/4), rel, max |dw|."""
    x = 0.25 * L
    w_ex = simply_supported_field(x)
    rows = []
    for n_el in (2, 4, 8, 16):
        u = simply_supported(n_el)
        w_fem = evaluate_w(u, x, L)
        rel = abs(w_fem - w_ex) / abs(w_ex)
        max_abs = 0.0
        for i in range(101):
            xi = L * i / 100.0
            dw = abs(evaluate_w(u, xi, L) - simply_supported_field(xi))
            if dw > max_abs:
                max_abs = dw
        rows.append((n_el, L / n_el, w_fem, w_ex, rel, max_abs))
    return rows


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    u = cantilever()
    w_tip = u[2]
    w_ref = cantilever_tip()
    print(f"cantilever  L = {L} m  E = {E:.3e} Pa  I = {I} m^4  P = {P} N")
    print(f"w_tip (1 element) = {w_tip:.12e} m")
    print(f"PL^3/3EI          = {w_ref:.12e} m")
    print(f"relative error    = {abs(w_tip - w_ref) / abs(w_ref):.6e}")

    u_ss = simply_supported()
    w_mid = evaluate_w(u_ss, 0.5 * L, L)
    w_mid_ref = simply_supported_midspan()
    print(f"simply-supported  q = {Q} N/m")
    print(f"w_mid (2 elements)     = {w_mid:.12e} m")
    print(f"5qL^4/384EI            = {w_mid_ref:.12e} m")
    print(f"relative error         = {abs(w_mid - w_mid_ref) / abs(w_mid_ref):.6e}")

    rows = mesh_sweep()
    print("n_el  h       w_fem(L/4)        w_exact(L/4)      rel          max|dw|     factor")
    last = None
    hs = []
    errs = []
    for n_el, h, w_fem, w_ex, rel, max_abs in rows:
        factor = (last / max_abs) if last is not None else float("nan")
        print(
            f"{n_el:4d}  {h:.4f}  {w_fem:.12e}  {w_ex:.12e}  {rel:.6e}  {max_abs:.6e}  {factor:.4f}"
        )
        last = max_abs
        hs.append(h)
        errs.append(max_abs)

    plot = Path(__file__).resolve().parent.parent / "error_vs_h.png"
    plt.loglog(hs, errs, "o-")
    plt.xlabel("h (m)")
    plt.ylabel("max |w_fem - w_exact| (m)")
    plt.title("SS uniform q, error vs h")
    plt.tight_layout()
    plt.savefig(plot)
    print(f"plot {plot}")
