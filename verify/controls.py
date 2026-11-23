"""Must-fail controls for the verifier.  Each must print PASS.
  C1 a NaN lower bound counts as a failing cell (_adaptive and _adaptive1d)
  C2 a non-finite or out-of-domain endpoint raises instead of being clipped
  C3 a perturbed node (G raised by 3 at node k) makes region 1 fail at small depth (run with --node)"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fei_cert23 as F
from fei_cert import Iv
import fei_cert_r1x as R

def say(name, ok): print("%s %s" % ("PASS" if ok else "FAIL", name), flush=True)

x0 = np.array([1.0, 2.0]); x1 = np.array([2.0, 3.0])
ok = F._adaptive(lambda X, Y: (np.array([1.0, np.nan])[:len(X.lo)] if len(X.lo) == 2 else np.full(len(X.lo), np.nan),
                               np.ones(len(X.lo)), np.ones(len(X.lo))), x0, x1, x0, x1, 2, "nan-control")
say("C1 NaN cell fails (2-D)", ok is False)
ok = F._adaptive1d(lambda A: np.where(A.lo > 1.5, np.nan, 1.0), x0, x1, 2, "nan-control-1d")
say("C1 NaN cell fails (1-D)", ok is False)

cert, C0 = R.load_cert(sys.argv[1], 0.01)
for bad, label in ((np.nan, "NaN"), (np.inf, "inf"), (0.5, "x=0.5 below the domain")):
    try:
        F.slope_range(cert, np.array([bad]), np.array([2.0])); raised = False
    except ValueError:
        raised = True
    say("C2 slope_range raises on %s" % label, raised)
    try:
        cert.Gnum_range(Iv(np.array([bad]), np.array([2.0]))); raised = False
    except ValueError:
        raised = True
    say("C2 Gnum_range raises on %s" % label, raised)

if "--node" in sys.argv:
    # C3, paired: the same region-1 cells, the same depth.  The true certificate must pass on them, the certificate
    # with G raised by 3 at node k must fail, and at a point of those cells an UPPER bound of Phi under the perturbed
    # G must be negative (so the failure is the verifier seeing a real violation, not a depth limit).
    d = np.load(sys.argv[1]); xs = d["xs"]; k = int(np.searchsorted(xs, 6.0))
    n = len(xs) - 1; I, J = np.triu_indices(n)
    s = (xs[I] + xs[J] <= xs[k]) & (xs[k] <= xs[I + 1] + xs[J + 1]); I, J = I[s][:8], J[s][:8]
    G = d["G"].copy(); G[k] += 3.0
    pc = F.Cert(xs, G, 0.01, tail_slope=(G[-1] - F.BSTAR) / float(xs[-1])); F.build_slope_tables(pc)
    args = (xs[I], xs[I + 1], xs[J], xs[J + 1], 12)
    ok_true = F._adaptive(lambda X, Y: R.phi1_lower_exact(cert, X, Y), *args, "C3 true")
    ok_pert = F._adaptive(lambda X, Y: R.phi1_lower_exact(pc, X, Y), *args, "C3 perturbed")
    x = np.array([xs[I[0]]]); y = np.array([xs[k] - xs[I[0]]])
    if not (xs[J[0]] <= y[0] <= xs[J[0] + 1]):
        y = np.array([xs[J[0]]]); x = np.array([xs[k] - xs[J[0]]])
    X, Y = Iv(x, x), Iv(y, y)
    M = Iv(0.01) * (1 - R.exp2neg(X)) * (1 - R.exp2neg(Y))
    up = (pc.Gnum_range(X) + pc.Gnum_range(Y) - pc.Gnum_range(X + Y) + M - R.E_iv(X, Y)).hi
    say("C3 true certificate passes the %d cells at depth 12" % len(I), ok_true is True)
    say("C3 node %d (x=%.4f) raised by 3: the same cells fail" % (k, xs[k]), ok_pert is False)
    say("C3 perturbed Phi at (%.6f, %.6f) has upper bound %.4f < 0" % (x[0], y[0], up[0]), bool(up[0] < 0))
