"""
closed_form.py

Hibbeler tip PL^3 / 3EI. Kept separate so k cannot check itself.
"""

# Same numbers I wrote down for the paper sitting. w points the same way as P.
L = 1.0  # beam length (m), wall to free end
E = 200.0e9  # Young's modulus (Pa). Steel-ish, not a specific alloy.
I = 1.0e-6  # second moment of area (m^4)
P = 1000.0  # tip load (N)
Q = 1000.0  # simply-supported uniform line load (N/m). Not P.


def cantilever_tip(P=P, length=L, E=E, I=I):
    """Wall at x = 0, tip load P. Returns PL^3 / 3EI."""
    return P * length**3 / (3.0 * E * I)


def simply_supported_midspan(q=Q, length=L, E=E, I=I):
    """Pins at both ends, uniform q. Returns 5 q L^4 / 384 EI."""
    return 5.0 * q * length**4 / (384.0 * E * I)


def simply_supported_field(x, q=Q, length=L, E=E, I=I):
    """Pins at both ends, uniform q. w(x) = (q x / 24 EI) (L^3 - 2 L x^2 + x^3)."""
    return (q * x / (24.0 * E * I)) * (length**3 - 2.0 * length * x**2 + x**3)


if __name__ == "__main__":
    w = cantilever_tip()
    print(f"PL^3 / (3EI) = {w:.12e} m")
    print(f"hand check   = {1.0 / 600.0:.12e} m")
    w_mid = simply_supported_midspan()
    print(f"5qL^4 / (384 EI) = {w_mid:.12e} m")
    print(f"field at L/2     = {simply_supported_field(0.5 * L):.12e} m")
