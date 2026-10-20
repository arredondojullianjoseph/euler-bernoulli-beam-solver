"""Hibbeler checks."""

from beam import cantilever, evaluate_w, simply_supported
from closed_form import L, cantilever_tip, simply_supported_midspan

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
