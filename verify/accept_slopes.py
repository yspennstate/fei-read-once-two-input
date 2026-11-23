"""Acceptance 1: every stored slope interval contains the exact secant of the dyadic nodes; report how many float
slopes miss it, and the first-piece numbers."""
import sys, os
from fractions import Fraction
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fei_cert import Cert

d = np.load(sys.argv[1]); xs, G = d["xs"], d["G"]; K = float(xs[-1])
BSTAR = 8 - 4 / np.log(2)
c = Cert(xs, G, 0.01, tail_slope=(G[-1] - BSTAR) / K)
miss_float = inside = widen = 0
for i in range(len(xs) - 1):
    s = (Fraction(float(G[i + 1])) - Fraction(float(G[i]))) / (Fraction(float(xs[i + 1])) - Fraction(float(xs[i])))
    assert Fraction(float(c.slopes_lo[i])) <= s <= Fraction(float(c.slopes_hi[i])), i
    inside += 1
    miss_float += Fraction(float(c.slopes[i])) != s
    widen += c.slopes_lo[i] < c.slopes_hi[i]
s0 = (Fraction(float(G[1])) - Fraction(float(G[0]))) / (Fraction(float(xs[1])) - Fraction(float(xs[0])))
print("pieces", len(xs) - 1, "intervals containing the exact secant", inside, "float slopes != exact", miss_float,
      "intervals of positive width", widen)
print("first piece exact", s0, "float", Fraction(float(c.slopes[0])), "interval", c.slopes_lo[0], c.slopes_hi[0])
neg = [i for i in range(len(xs) - 1) if c.slopes_hi[i] < 0]
print("negative pieces", neg, "min slope on [K-1,inf)", min(c.slopes_lo[xs[1:] > K - 1].min(), c.tail_slope), "tail", c.tail_slope)
