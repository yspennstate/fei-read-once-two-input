# The Fourier entropy-influence constant of read-once formulas over two-input gates

Yitzchak Shmalo

For every read-once formula f over gates of fan-in two, with negations, on independent uniform bits,
H(f) <= 6.4573 I(f). Finite truncations of Hod's construction show that no constant below 6.4547837 is possible for
this class.

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `verify/`: the certificate `cert_vi_h24k40.npz` and the programs that check it.
- `logs/`: the output of one run of those programs.

## Running the verification

Python 3 with numpy. From `verify/`:

    sha256sum -c SHA256SUMS.txt
    python3 run_cert_all_x.py cert_vi_h24k40.npz 0.01 30 run
    python3 controls.py cert_vi_h24k40.npz --node

`run_cert_all_x.py` runs the six parts of the certificate, each as its own process, and prints
`ALL SIX PARTS VERIFIED` only if every part prints its results line with no `False` in it and exits with code 0.
The parts are the tail lemmas and constants (`lemmas`), the one-dimensional checks T2 and T3, and the
two-dimensional regions R1, R2 and R3. All bounds are computed in interval arithmetic with outward rounding.
`controls.py` holds the must-fail controls; each prints PASS. The run in `logs/` took about three minutes on one core.

Other programs:

| file | what it checks |
|---|---|
| `accept_slopes.py` | every stored slope interval contains the exact secant of the dyadic nodes (rational arithmetic) |
| `exact_const.py` | the exact rational value of the certified constant against 64573/10000 |
| `xor_premises.py` | the premises of the proof of the XOR lemma on the stored nodes |
| `xor_lemma_check.py` | the XOR lemma on a grid, in floating point (a numerical check, not a proof) |
| `ctl_straddle.py` | a control for the range enclosure on pieces where the slope is negative |
| `xor_break.py`, `xor_climb.py` | random and hill-climbing search over small read-once formulas for a ratio H/I above 6.4573 |

`fei_cert.py`, `fei_cert23.py`, `fei_cert_r1x.py` and `run_r2x.py` are the verifier itself and are called by
`run_cert_all_x.py`.

## The run in `logs/`

| file | contents |
|---|---|
| `run_all.log` | hash check of every program and the certificate, one line per part, the summary line |
| `certx_run_lemmas.out` | the tail lemmas and the constants |
| `certx_run_T2.out`, `certx_run_T3.out` | the one-dimensional checks T2 and T3 |
| `certx_run_R1x.out`, `certx_run_R2x.out`, `certx_run_R3.out` | the two-dimensional checks R1, R2, R3 |
| `certx_run_SUMMARY.json` | exit code and time of each part |
| `controls.log` | the must-fail controls |
| `timing.txt` | time and memory of the run |
| `SHA256SUMS.txt` | sha256 of the files above |

Results, as printed:

- all six parts verified, exit code 0;
- R1: 1,620,900 initial cells, 1,699,336 evaluated, deepest level 6;
- R2: 3,240,000 initial cells, 15,825,700 evaluated, deepest level 24;
- R3: 1,656,900 initial cells, 23,299,826 evaluated, deepest level 28;
- T2: deepest level 5, smallest lower bound 2.588e-03;
- T3: deepest level 2, smallest lower bound 1.739e-03;
- T1 margin 1.000e-02; region 2 for v < 2^-K, lower bound 120.486; region 3 for u < 2^-60, largest derivative
  198.99, lower bound 7.98; region 3 for u, v <= 2^-K, smallest slope on [K-1, inf) 3.146618, lower bound 4.539611;
  |G'| <= 30.085, |Psi'| <= 21.892;
- controls: 11 PASS.
