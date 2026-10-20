"""
beam.py

Cantilever: one element vs Hibbeler PL^3 / 3EI.
Simply-supported: two equal elements so L/2 is a node vs 5qL^4 / 384EI.
"""

import numpy as np

from closed_form import E, I, L, P, Q, cantilever_tip, simply_supported_midspan
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
    """Deflection at x. 4-DOF: one element. 6-DOF: two equal elements on [0, length]."""
    if len(u) == 4:
        return float(hermite_shape(x / length, length) @ u)
    le = 0.5 * length
    if x <= le:
        return float(hermite_shape(x / le, le) @ u[:4])
    return float(hermite_shape((x - le) / le, le) @ u[2:6])


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


def simply_supported():
    """Pins on w at both ends. Uniform q. Two equal elements; midspan is nodal w."""
    n_el = 2
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


if __name__ == "__main__":
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
