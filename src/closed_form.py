"""
closed_form.py

Hibbeler tip PL^3 / 3EI. Kept separate so k cannot check itself.
"""

# Same numbers I wrote down for the paper sitting. w points the same way as P.
L = 1.0  # beam length (m), wall to free end
E = 200.0e9  # Young's modulus (Pa). Steel-ish, not a specific alloy.
I = 1.0e-6  # second moment of area (m^4)
P = 1000.0  # tip load (N)


def cantilever_tip(P=P, length=L, E=E, I=I):
    """Wall at x = 0, tip load P. Returns PL^3 / 3EI."""
    return P * length**3 / (3.0 * E * I)


if __name__ == "__main__":
    w = cantilever_tip()
    print(f"PL^3 / (3EI) = {w:.12e} m")
    print(f"hand check   = {1.0 / 600.0:.12e} m")
