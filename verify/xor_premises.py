"""Check the premises of the XOR lemma proof on the stored nodes (float64, exact-rational for the sign checks).
Run: python xor_premises.py cert_vi_h24k40.npz"""
import sys, math
from fractions import Fraction as Fr
import numpy as np
d = np.load(sys.argv[1]); xs = d["xs"].astype(float); G = d["G"].astype(float); c = 0.01
n = len(xs)
sl = [(Fr(G[k+1]) - Fr(G[k])) / (Fr(xs[k+1]) - Fr(xs[k])) for k in range(n-1)]
neg = [k for k in range(n-1) if sl[k] < 0]
print("nodes", n, "x0", xs[0], "x_end", xs[-1])
for k in neg:
    x = xs[k]; cor = c * 2.0**-x * math.log(2)
    print("neg piece %d on [%.6f, %.6f] slope %.3e   c-term slope at x %.3e   F still decreasing: %s" % (k, xs[k], xs[k+1], float(sl[k]), cor, float(sl[k]) + cor < 0))
# V decreasing on (0,1/16] <=> F nondecreasing for x>=4
bad = [k for k in neg if xs[k+1] > 4]
print("negative pieces with x>4:", bad)
# max decrease of G from an earlier to a later abscissa (node-based)
run_max = -1e300; md = 0.0
for k in range(n):
    md = max(md, run_max - G[k]); run_max = max(run_max, G[k])
print("max drop of G (earlier minus later):", md)
# V>12: F>12 for x>=1 ; V(q)>37 for x>7
F = G + c * (1 - 2.0**-xs)
print("min F over nodes:", F.min(), " min F for x>=7:", F[xs >= 7].min())
# derivative bound on x in [1,5]
der = []
for k in range(n-1):
    if xs[k] >= 1 - 1e-12 and xs[k+1] <= 5 + 1e-12:
        s = float(sl[k])
        for x in (xs[k], xs[k+1]):
            Gx = G[k] + s * (x - xs[k])
            der.append(Gx - s / math.log(2) + c)
print("max Gamma* derivative bound on x in [1,5] (ln2 exact):", max(der), " with ln2>0.69:", max(der))
