"""Float check of the XOR lemma for the read-once certificate potential:
    Gamma*(p+q-2pq) <= Gamma*(p) + Gamma*(q)  for p, q in [0,1],
Gamma*(p) = p F(log2(1/p)) for p <= 1/2, F(x) = G(x) + c(1-2^-x), G piecewise linear through the npz nodes, linear tail for
x > 40 with slope (G(40) - 8 + 4/ln2)/40; reflected across 1/2, zero at 0 and 1.  Float64 only: a numerical check, not a proof.
Run: python xor_lemma_check.py cert_vi_h24k40.npz"""
import sys, hashlib, math
import numpy as np
fn = sys.argv[1]
print("npz sha256", hashlib.sha256(open(fn, "rb").read()).hexdigest())
d = np.load(fn); xs = d["xs"].astype(float); G = d["G"].astype(float); c = 0.01
sig = (G[-1] - 8 + 4 / math.log(2)) / xs[-1]
print("nodes", len(xs), "x range", xs[0], xs[-1], "tail slope", sig, "C", d["C"], "K", d["K"], "M", d["M"])

def Gx(x):
    x = np.asarray(x, float)
    return np.where(x <= xs[-1], np.interp(x, xs, G), G[-1] + sig * (x - xs[-1]))

def gam(p):
    p = np.asarray(p, float)
    s = np.minimum(p, 1 - p)
    out = np.zeros_like(s)
    m = s > 0
    x = -np.log2(s[m])
    out[m] = s[m] * (Gx(x) + c * (1 - 2.0 ** -x))
    return out

print("Gamma*(1/2) =", float(gam(np.array([0.5]))[0]), "(leaf needs < 6.4573)")
# grid: x = log2(1/p) from 1 to 80 dense near 1, plus linear grid in p near 1/2 and beyond 1/2
x1 = np.concatenate([np.linspace(1, 6, 2001), np.linspace(6, 80, 1500)])
P = np.unique(np.concatenate([2.0 ** -x1, np.linspace(1e-4, 1 - 1e-4, 3001)]))
worst = (np.inf, None, None); worst_rel = (np.inf, None, None); neg = 0; tot = 0
for i in range(0, len(P), 400):
    p = P[i:i + 400, None]; q = P[None, :]
    r = p + q - 2 * p * q
    D = gam(p) + gam(q) - gam(r)
    tot += D.size; neg += int((D < -1e-12).sum())
    k = np.unravel_index(np.argmin(D), D.shape)
    if D[k] < worst[0]: worst = (float(D[k]), float(p[k[0], 0]), float(q[0, k[1]]))
    rel = D / np.minimum(np.minimum(p, 1 - p), np.minimum(q, 1 - q))
    k = np.unravel_index(np.argmin(rel), rel.shape)
    if rel[k] < worst_rel[0]: worst_rel = (float(rel[k]), float(p[k[0], 0]), float(q[0, k[1]]))
print("pairs", tot, "negative (< -1e-12)", neg)
print("min D = %.6g at p=%.6g q=%.6g" % worst)
print("min D/min-side = %.6g at p=%.6g q=%.6g" % worst_rel)
print("PASS" if neg == 0 else "FAIL")
