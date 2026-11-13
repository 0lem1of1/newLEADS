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

### 2. Model 0 envelope derivative lacked `sign(eta)`

**Where:** `leads_laser.cpp:46-47`, and its copy in `cuda/leads_boris.cuh:79-80`.
Only Model 0 (plane wave) uses this term.

**Problem.** The field is meant to be `E = -dA/dt` with
`A = a0 * g(eta) * cos(eta)` (and `sin` for the other component), where
`eta = t - z + phase` and `g = exp(-alpha*|eta|^n)`. The code wrote the
envelope term as

```cpp
derivative = alpha * order * pow( abs(eta),order-1) ;   // = -g'/g only for eta > 0
```

but `d/deta |eta|^n = n*|eta|^(n-1)*sign(eta)`. For `eta < 0` (the leading half
of the pulse) the term had the wrong sign. E and B shared the error, so `B = z x E`
and the invariant `gamma - p_z` still held, which hid it. The field was no
longer the derivative of a bounded potential, so `int E dt != 0` and the pulse
left a net transverse kick. Theory says a plane-wave pulse leaves zero net
transverse momentum.

The CUDA port copied the formula exactly, so CPU/GPU agreement could not
detect it, and neither could self-convergence.

**Fix.**

```cpp
derivative = alpha * order * pow( abs(eta),order-1) * copysign(1.0,eta) ;
```

The same change was made in `laser_profile<0>`, so LEADS and the harness still
match. When `PulseProfile = 0` (infinite plane wave), `alpha = 0` and the term
is zero either way.

**How it was checked** (scratch builds of the tree before and after the fix,
2026-11-13).

- *Against theory: net impulse.* 20 electrons, gamma = 500, head-on into an
  a0 = 10 pulse, with `z0 = 2.5*Tau_FWHM`, so the envelope is about 1e-30 at
  both ends of the window. The table gives the residual `max|p_x(final) -
  p_x(0)|/a0`:

  | Pulse (cycles) | Linear before | Linear after | Circular before | Circular after |
  |---|---|---|---|---|
  | 1 | 0.46 | 2.2e-7 | 0.33 | 5.6e-8 |
  | 2 | 0.079 | 1.2e-8 | 0.056 | 2.8e-9 |
  | 4 | 0.018 | 2.2e-14 | 0.013 | 1.2e-14 |
  | 10 | 0.0028 | 3.5e-14 | 0.0020 | 3.7e-14 |

  The 1-cycle leftover halves each time `dTau` is halved (2.2e-7, 1.1e-7,
  5.4e-8, 2.7e-8 at `--refine` 1, 2, 4, 8). It is the pusher's first-order time
  error (items 3 and 4), not the field. After the fix, `gamma` returns to its
  initial value to 1e-10, and before the fix it was off by up to 2e-5.
  `gamma - p_z` is conserved to 1e-10 or better in both builds.
- *Against the harness.* Reference case (5000 particles, 4 cycles, a0 = 10):
  Armadillo vs `bench_cpu` is `BIT-IDENTICAL`. For `bench_cpu` vs `bench_gpu`,
  the L2 gap is posi 2.6e-13 and momt 2.1e-10, the same as before the fix (2.6e-13
  and 2.4e-10). The worst relative error rose (velo 1e-8, was 7e-13) because
  the final transverse velocity is now close to zero. Absolute errors did not
  change.
- *Model 2 is untouched.* LEADS before vs after, Model 2 with `pola` 0 and 1, is
  `BIT-IDENTICAL`.
- *Convergence unchanged.* `convtest` order is 1.00 on posi/velo/momt before and
  after. That is the known first-order result (items 3 and 4).
- `validation/thomson_compton_check.py` was not rerun. It uses a0 = 0.01 and
  10 cycles, where the old residual kick is about 3e-3 of a0.

**After pulling this change.** Every Model 0 result changes. In the reference
case, Armadillo before vs after differs by up to 0.18 in `momt` and 7.6e-3 in
`posi`. Regenerate `TestCase/bench/state_final_arma.bin` and any Model 0 numbers
quoted in the reports. Model 2 data is still valid.

## Known, not yet fixed

| # | Problem | Where |
|---|---|---|
| 3 | Boris computes gamma from the momentum before the half E-kick instead of from `UM` | `leads_traj.cpp:20-23`, `cuda/leads_boris.cuh:167-170` |
| 4 | `--halfstep` resync kicks by `dTau/2` instead of `dTau` | `cuda/bench_cpu.cpp:179` |
