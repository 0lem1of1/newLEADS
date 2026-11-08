# Fixes to the original LEADS code

Line numbers refer to the tree at the commit that adds this file.

## Fixed

### 1. The time axis `TR` did not advance by `dTau`

**Where:** `leads_dyna.cpp:39-40` (the fix). Related code:

| Line | What it does |
|---|---|
| `leads_dyna.cpp:11` | `dTau = in.TimeStep`, the step the pusher takes |
| `leads_dyna.cpp:22` | `nT = int(2*t_shift/dTau) + 1` |
| `leads_dyna.cpp:81` | `CalForceMovePart(TR(i), j)`, so `TR` is the lab time handed to the laser |
| `leads_traj.cpp:152` | `rp = rp + beta*dTau`, so the particle moves by `dTau` per step (each radiation-reaction branch has its own copy at lines 51, 95, 118, 140) |
| `leads_radi.cpp:92` | `vT(np,i) = TR(i) - n.r`, the retarded time that goes into the spectrum phase |
| `leads_init.cpp:82-85` | `t_shift = z0/beta` for an electron beam (`GamaParticle > 1`) |

**Problem.** The pusher moves particles by `dTau` per step, but `TR` was built as

```cpp
TR = linspace(-t_shift, t_shift, nT);   // step = 2*t_shift/(nT-1)
```

so its step was `2*t_shift/(nT-1)`, not `dTau`. For an electron beam
`t_shift = z0/beta`, and `nT` truncates `2*t_shift/dTau` to an integer, so the
step comes out as about `dTau/beta`. The relative gap is `1/beta - 1` (about
`1-beta`) when `2*z0/dTau` is a whole number, and the truncation of `nT` can add
up to `dTau/(2*t_shift)` on top of that.

That gap is tiny for the trajectory (the laser phase the particle sees is off by
about `(1-beta)/2`), but the radiation uses the retarded time `t - n.r`. For
emission along the electron's direction one step advances it by
`dTau*(1-beta)` in reality, and by `TR_step - beta*dTau` in the code. That
difference is comparable to `1-beta` itself, so the retarded time is stretched
by a large factor and the spectrum peak lands far below the true Thomson
frequency.

CPU/GPU agreement and self-convergence cannot see it: the gap does not depend
on `dTau`, and `bench_cpu`/`bench_gpu` already used the consistent grid.

**Fix.** Build the axis from `dTau` directly:

```cpp
TR = zeros(nT);
for(int i = 0; i < nT; i++) TR(i) = -t_shift + i*dTau;   // leads_dyna.cpp:39-40
```

This is the same arithmetic as `cuda/bench_cpu.cpp:109` and
`cuda/bench_gpu.cu:161`. The step is now exactly `TimeStep`. The cost is that
the window ends at `-t_shift + (nT-1)*dTau`, which can be up to one `dTau` short
of `+t_shift`.

**How it was checked.**

- *Against the harness.* LEADS (500 particles, Model 0) replayed through
  `bench_cpu` and compared with `statediff`. Before: L2 gap posi 6.8e-8, momt
  1.8e-5. After: `BIT-IDENTICAL` on all five arrays.
- *Against theory.* `validation/thomson_compton_check.py` runs one on-axis
  electron head-on into an a0 = 0.01, 10-cycle pulse and compares the spectrum
  peak at theta = 180 deg with the Thomson frequency `(1+beta)/(1-beta)`:

  | gamma | Peak error before | Peak error after |
  |---|---|---|
  | 50 | -13.27% | -0.003% |
  | 100 | -44.26% | -0.003% |
  | 500 | -94.03% | -0.003% |

  The "before" error depends on the truncation in `nT`; with `z0 = 20` and
  `dTau = 5e-3` (where `2*z0/dTau` is a whole number) the same arithmetic gives
  0.5x theory for gamma = 100 and 500 and 0.73x for gamma = 50. The leftover
  -0.003% is the same for every gamma and pulse length, and is about the size of
  the a0^2/2 redshift. Halving `dTau` changes the spectrum by at most 1.4e-5 in
  `ln(y)`.
- LEADS has no electron recoil, so Thomson is the correct reference. The
  Compton edge (printed alongside) lies 0.061%, 0.121% and 0.603% below it at
  gamma = 50, 100 and 500.

**After pulling this change.** Rebuild with `make all`. Reference data generated
before it (including `TestCase/bench/state_final_arma.bin` and the L2 gaps
against Armadillo quoted in the reports) used the old grid and must be
regenerated. Those gaps should drop to zero.

## Known, not yet fixed

| # | Problem | Where |
|---|---|---|
| 2 | Model 0 envelope derivative lacks `sign(eta)`, so E != -dA/dt and the pulse leaves a net impulse | `leads_laser.cpp:46`, copied into `cuda/leads_boris.cuh` |
| 3 | Boris computes gamma from the momentum before the half E-kick instead of from `UM` | `leads_traj.cpp:20-23`, `cuda/leads_boris.cuh:167-170` |
| 4 | `--halfstep` resync kicks by `dTau/2` instead of `dTau` | `cuda/bench_cpu.cpp:179` |
