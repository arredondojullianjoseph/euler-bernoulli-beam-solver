"""
beam.py

One-element cantilever vs Hibbeler PL^3 / 3EI.
"""

import numpy as np

from closed_form import E, I, L, P, cantilever_tip
from stiffness import element_stiffness


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


if __name__ == "__main__":
    u = cantilever()
    w_tip = u[2]
    w_ref = cantilever_tip()
    print(f"cantilever  L = {L} m  E = {E:.3e} Pa  I = {I} m^4  P = {P} N")
    print(f"w_tip (1 element) = {w_tip:.12e} m")
    print(f"PL^3/3EI          = {w_ref:.12e} m")
    print(f"relative error    = {abs(w_tip - w_ref) / abs(w_ref):.6e}")
