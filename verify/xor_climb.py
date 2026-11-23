"""Hill-climb version of xor_break.py: mutate one random subtree (re-sample it, or flip a gate type / negation), keep if H/I rises.
Several restarts, including from pure-AND/OR seeds so XOR insertions are tried against the AND/OR extremals.
Run: python xor_climb.py SEED NMAX STEPS RESTARTS"""
import sys, random
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], "0"] + sys.argv[3:]
import importlib.util, os
spec = importlib.util.spec_from_file_location("xb", os.path.join(os.path.dirname(os.path.abspath(__file__)), "xor_break.py"))
xb = importlib.util.module_from_spec(spec); spec.loader.exec_module(xb)
STEPS, RS = int(sys.argv[4]), int(sys.argv[5]); NMAX = int(sys.argv[2]); rnd = xb.rnd
def nodes(t, path=()):
    yield path
    if t[0] != "x":
        yield from nodes(t[2], path + (2,)); yield from nodes(t[3], path + (3,))
def get(t, p):
    for i in p: t = t[i]
    return t
def put(t, p, s):
    if not p: return s
    l = list(t); l[p[0]] = put(t[p[0]], p[1:], s); return tuple(l)
def mutate(t):
    p = rnd.choice(list(nodes(t))); s = get(t, p); m = rnd.random()
    if s[0] != "x" and m < .3: s = (rnd.choice("AOX"), s[1], s[2], s[3])
    elif s[0] != "x" and m < .45: s = (s[0], not s[1], s[2], s[3])
    else:
        room = NMAX - xb.leaves(t) + xb.leaves(s)
        s = xb.rand_tree(rnd.randint(1, max(1, room)))
    return put(t, p, s)
gbest, gt = 0, None
for r in range(RS):
    t = xb.rand_tree(rnd.randint(2, NMAX)); cur = xb.ratio(t)
    for _ in range(STEPS):
        u = mutate(t)
        if xb.leaves(u) > NMAX: continue
        v = xb.ratio(u)
        if v >= cur: t, cur = u, v
    has_x = "'X'" in repr(t)
    print("restart %d: %.6f  XOR in best: %s  leaves %d" % (r, cur, has_x, xb.leaves(t)))
    if cur > gbest: gbest, gt = cur, t
print("GLOBAL max H/I %.6f (bound 6.4573)" % gbest); print(gt)
