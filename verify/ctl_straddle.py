"""Control C4: Gnum_range on intervals straddling the negative-slope pieces 510, 511, 767, 895 contains every
sampled value of the piecewise-linear G_num (exact rational evaluation at 2001 points) and its lower end is at or below
the true minimum (an interior node).  Run: python ctl_straddle.py cert_vi_h24k40.npz"""
import sys
from fractions import Fraction as Fr
import numpy as np
import fei_cert as fc

import fei_cert23 as F
d = np.load(sys.argv[1])
cert = F.Cert(d["xs"], d["G"], 0.01, tail_slope=(d["G"][-1] - F.BSTAR) / float(d["xs"][-1])); F.build_slope_tables(cert)
if cert is None:
    raise SystemExit("no Cert class in fei_cert")
xs, G = cert.xs, cert.G
ok = True
for k in (510, 511, 767, 895):
    a, b = float(xs[k - 1]) + 1e-9, float(xs[k + 2]) - 1e-9      # covers pieces k-1, k, k+1 with interior nodes
    r = cert.Gnum_range(fc.Iv(np.array([a]), np.array([b])))
    lo, hi = float(r.lo[0]), float(r.hi[0])
    vals = []
    for t in np.linspace(a, b, 2001):
        j = int(np.searchsorted(xs, t, side="right") - 1)
        v = Fr(float(G[j])) + (Fr(float(t)) - Fr(float(xs[j]))) * (Fr(float(G[j + 1])) - Fr(float(G[j]))) / (Fr(float(xs[j + 1])) - Fr(float(xs[j])))
        vals.append(v)
    vmin, vmax = min(vals), max(vals)
    good = Fr(lo) <= vmin and vmax <= Fr(hi)
    ok &= good
    print("piece %d slope %+.4f  range [%.9f, %.9f]  sampled [%.9f, %.9f]  %s" % (
        k, (float(G[k + 1]) - float(G[k])) / (float(xs[k + 1]) - float(xs[k])), lo, hi, float(vmin), float(vmax), "PASS" if good else "FAIL"))
print("C4", "PASS" if ok else "FAIL")
