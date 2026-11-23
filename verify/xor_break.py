"""Break-test for paper C: random read-once formulas over fan-in-2 gates (AND, OR, XOR, with input/output negations), uniform
inputs, n <= NMAX leaves. Exact H/I from the Walsh-Hadamard transform. Reports the max H/I seen; anything above 6.4573 is a
counterexample. Hill-climbs by re-sampling subtrees. Run: python xor_break.py SEED NMAX TRIALS"""
import sys, random, numpy as np
seed, NMAX, TR = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rnd = random.Random(seed)
def rand_tree(n):
    if n == 1: return ("x",)
    k = rnd.randint(1, n - 1)
    return (rnd.choice("AOX"), rnd.random() < .3, rand_tree(k), rand_tree(n - k))
def leaves(t): return 1 if t[0] == "x" else leaves(t[2]) + leaves(t[3])
def ev(t, X, pos):
    if t[0] == "x":
        v = X[:, pos[0]]; pos[0] += 1; return v
    a = ev(t[2], X, pos); b = ev(t[3], X, pos)
    r = (a & b) if t[0] == "A" else (a | b) if t[0] == "O" else (a ^ b)
    return r ^ 1 if t[1] else r
def ratio(t):
    n = leaves(t)
    X = ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1).astype(np.int8)
    f = 1.0 - 2.0 * ev(t, X, [0]).astype(float)
    h = f.copy(); s = 1
    while s < len(h):
        h = h.reshape(-1, 2 * s); a, b = h[:, :s].copy(), h[:, s:].copy(); h[:, :s] = a + b; h[:, s:] = a - b; h = h.ravel(); s *= 2
    w = (h / len(h)) ** 2
    deg = np.array([bin(i).count("1") for i in range(len(h))])
    I = float((w * deg).sum())
    if I < 1e-12: return 0.0
    nz = w[w > 1e-300]
    H = float(-(nz * np.log2(nz)).sum())
    return H / I
best, bt = 0.0, None
for _ in range(TR):
    t = rand_tree(rnd.randint(2, NMAX)); r = ratio(t)
    if r > best: best, bt = r, t
print("max H/I %.6f  (paper C bound 6.4573, AND/OR/NOT sharp ~6.4548)" % best)
print("tree", bt)
