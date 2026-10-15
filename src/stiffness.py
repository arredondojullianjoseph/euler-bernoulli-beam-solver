"""
stiffness.py

One-element cubic Hermite k from k_ij = integral EI N_i'' N_j'' dx.
"""

import numpy as np

from closed_form import E, I, L


def second_deriv(xi, length):
    """d²N/dx² at xi = x/length. Order: w_left, theta_left, w_right, theta_right."""
    # N1 = 1 - 3 xi^2 + 2 xi^3
    #   d²/dξ² = -6 + 12 xi, then / length^2
    N1pp = (-6.0 + 12.0 * xi) / length**2
    # N2 = length * (xi - 2 xi^2 + xi^3)
    #   d²/dξ² = length * (-4 + 6 xi), then / length^2 → / length
    N2pp = (-4.0 + 6.0 * xi) / length
    # N3 = 3 xi^2 - 2 xi^3
    #   d²/dξ² = 6 - 12 xi, then / length^2
    N3pp = (6.0 - 12.0 * xi) / length**2
    # N4 = length * (-xi^2 + xi^3)
    #   d²/dξ² = length * (-2 + 6 xi), then / length^2 → / length
    N4pp = (-2.0 + 6.0 * xi) / length
    return [N1pp, N2pp, N3pp, N4pp]


def element_stiffness(EI, length):
    """4x4 k. Two Gauss points; exact because N'' is linear."""
    # Gauss-Legendre on [0, 1] (mapped from ±1/sqrt(3) on [-1, 1]).
    points = [0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)]
    weight = 0.5

    k = np.zeros((4, 4))
    for xi in points:
        Npp = second_deriv(xi, length)
        for i in range(4):
            for j in range(4):
                # EI * length * w * N_i'' * N_j'' at this sample point
                k[i, j] += EI * length * weight * Npp[i] * Npp[j]
    return k


if __name__ == "__main__":
    EI = E * I
    k = element_stiffness(EI, L)
    print(k)
    print("k[0,0] / (EI/L^3) should be 12:")
    print(k[0, 0] / (EI / L**3))
