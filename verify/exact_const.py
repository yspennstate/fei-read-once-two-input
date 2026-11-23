"""Exact rational value of C_cert = G[0]/2 + c/4 for the stored node and binary64 c, against 64573/10000."""
import sys, math
from fractions import Fraction
import numpy as np
d = np.load(sys.argv[1]); G0 = Fraction(float(d["G"][0])); T = Fraction(64573, 10000)
for c in (0.01, math.nextafter(0.01, 0)):
    C = G0 / 2 + Fraction(c) / 4
    print("c=%r  C_cert - 6.4573 = %s (%.3e)  -> %s" % (c, C - T, float(C - T), "<= 6.4573" if C <= T else "ABOVE 6.4573"))
print("G0/2 =", float(G0 / 2), " largest c allowed =", float((T - G0 / 2) * 4))
