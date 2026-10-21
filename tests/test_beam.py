"""Hibbeler checks and mesh-order at L/4."""

from beam import cantilever, evaluate_w, mesh_sweep, simply_supported
from closed_form import L, cantilever_tip, simply_supported_field, simply_supported_midspan

REL_TOL = 1e-12


def test_cantilever_tip():
    u = cantilever()
    w_tip = u[2]
    w_ref = cantilever_tip()
    err = abs(w_tip - w_ref) / abs(w_ref)
    print(f"cantilever w_tip = {w_tip:.16e}  Hibbeler = {w_ref:.16e}  rel = {err:.6e}")
    assert err < REL_TOL


def test_simply_supported_midspan():
    u = simply_supported()
    w_mid = evaluate_w(u, 0.5 * L, L)
    w_ref = simply_supported_midspan()
    err = abs(w_mid - w_ref) / abs(w_ref)
    print(f"ss w_mid = {w_mid:.16e}  Hibbeler = {w_ref:.16e}  rel = {err:.6e}")
    assert err < REL_TOL


def test_mesh_order_quarter_span():
    field_mid = simply_supported_field(0.5 * L)
    hibbeler_mid = simply_supported_midspan()
    assert abs(field_mid - hibbeler_mid) / abs(hibbeler_mid) < REL_TOL

    rows = mesh_sweep()
    print("n_el  h       w_fem(L/4)        w_exact(L/4)      rel          max|dw|     factor")
    last = None
    maxes = []
    for n_el, h, w_fem, w_ex, rel, max_abs in rows:
        factor = (last / max_abs) if last is not None else float("nan")
        print(
            f"{n_el:4d}  {h:.4f}  {w_fem:.16e}  {w_ex:.16e}  {rel:.6e}  {max_abs:.6e}  {factor:.4f}"
        )
        last = max_abs
        maxes.append(max_abs)

    # 2-element L/4 is ~1.8%, not machine eps. L/4 is a node after that.
    assert 1e-2 < rows[0][4] < 3e-2
    assert all(m > 0.0 for m in maxes)
    for i in range(len(maxes) - 1):
        assert maxes[i + 1] < maxes[i]
