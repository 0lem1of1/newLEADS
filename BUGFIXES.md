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
  error (item 3, fixed below), not the field. After the fix, `gamma` returns to its
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

### 3. The Boris rotation used gamma from before the half E-kick

**Where:** `leads_traj.cpp:19-31` and `cuda/leads_boris.cuh` `boris_step`.

**Problem.** The magnetic rotation is scaled by `1/gamma`, which should be the
mid-step value: `gamma = sqrt(1+|UM|^2)` with `UM = P + del_s*E` (Boris 1970).
The rotation keeps `|UM|` fixed, so this is also the value after the rotation.
LEADS took `gamma` from `P`, the momentum at the start of the step. That is off
by about `(q*dTau/2)(u.E)/gamma`, so every rotation had an O(dTau) relative
error. The result was a first-order scheme. It only matters where E and B are
both nonzero. Radiation-reaction terms use the separate `gamaOld`, which is
unchanged.

**Fix.** Compute `UM` first, then `gamaM = sqrt(1 + cdot(UM,UM))`, and use it
for `del`. `boris_step` makes the same change with the same operation order as
`cdot`, so Armadillo and the harness still match.

**How it was checked** (scratch builds before and after, 2026-11-15; "before"
is the tree with fixes 1 and 2).

- *Net impulse.* Same test as item 2 (20 electrons, gamma = 500, a0 = 10,
  z0 = 2.5*Tau_FWHM). `max|dp_x|/a0` is the sideways momentum left over; it
  should be zero:

  | Case | Before, refine 1, 2, 4, 8 | After, refine 1, 2, 4, 8 |
  |---|---|---|
  | Linear, 1 cycle | 2.2e-7, 1.1e-7, 5.4e-8, 2.7e-8 | 4e-15, 6e-15, 6e-15, 1e-14 |
  | Linear, 2 cycles | 1.2e-8, 5.7e-9, 2.9e-9, 1.4e-9 | 1e-14, 1e-14, 2e-14, 2e-14 |
  | Circular, 1 cycle | 5.6e-8, 2.8e-8, 1.4e-8, 7.0e-9 | 3e-15, 4e-15, 7e-15, 1e-14 |
  | Circular, 2 cycles | 2.8e-9, 1.4e-9, 6.9e-10, 3.5e-10 | 7e-15, 8e-15, 9e-15, 2e-14 |

  The leftover halves with each refinement before the fix (first order) and
  is at round-off after it. So this bug was the whole of the item-2 leftover.
  The 4-cycle cases were already at round-off.
- *Convergence (`convtest --order 2 --tol 0.25`).* Without `--halfstep`, the
  order on posi/velo/momt is still 1.00 after the fix, for Model 0 and for
  Model 2 with `pola` 0 and 1. The final momentum is staggered half a step
  behind the position, and that offset differs on every grid. That is item 4,
  and the fix to item 3 cannot remove it. To separate the two, item 4's
  `bench_cpu.cpp:179` kick was changed from `dTau/2.0` to `dTau` in a scratch
  copy only (the repo is unchanged), run with `--halfstep`:

  | Case | Before: posi, velo, momt | After: posi, velo, momt |
  |---|---|---|
  | Model 0, reference | 2.03, 1.45, 1.06 | 2.00, 2.00, 2.00 (PASS) |
  | Model 2, `pola` 0 | 1.99, 1.64, 1.29 | 2.00, 2.00, 2.00 (PASS) |
  | Model 2, `pola` 1 | not run | 2.00, 2.00, 2.00 (PASS) |

  Neither fix alone gives order 2 on velo and momt. `egama` is not gated, and
  its error is now near round-off (5e-10), so its order is not meaningful.
- *Against the harness.* Armadillo vs `bench_cpu` is `BIT-IDENTICAL` for the
  Model 0 reference case. For Model 2 the L2 gap is posi 4e-14 and momt 3e-13
  (`pola` 0), as it was before. `bench_cpu` vs `bench_gpu` has L2 gaps of about
  2.6e-13 (posi) and 2.2e-10 (momt) for Model 0, and 1e-13 and 9e-11 for Model 2.
  Both are the same as before.
- *Plane-wave invariant.* `gamma - p_z` is conserved to 3e-14 (relative) in the
  Model 0 reference case.
- *Radiation reaction.* LEADS with `RadiationReaction = true` and `RRModel`
  0-3 (20 particles, 2-cycle pulse) runs and gives finite output. Momentum
  differs from before by up to 1.2e-7, which is about the O(dTau) size
  expected. The harness has no radiation reaction, so this only tests LEADS.

**After pulling this change.** Model 0 and Model 2 results both change slightly.
In the reference case Armadillo before vs after differs by up to 8.4e-8 in
`momt` and 3.4e-9 in `posi`, and in Model 2 by up to 1.0e-7 and 5.7e-9. Regenerate
`TestCase/bench/state_final_arma.bin` and the reference data for the reports.
Quoted CPU/GPU and Armadillo gaps will change at the 1e-7 level in `momt`.

### 4. The `--halfstep` resync kicked by `dTau/2` instead of `dTau`

**Where:** `cuda/bench_cpu.cpp`, the end-of-loop block under `if(halfstep)`.
Only the `--halfstep` path of `bench_cpu` is affected; LEADS, `bench_gpu` and the
default `bench_cpu` path are untouched.

**Problem.** `--halfstep` kicks momentum back to t0 - dTau/2 before the loop, so
after the loop it is p_{n-1/2}, half a step behind position. To put it back on
`t_final` the code took one more forward kick and averaged. It used
`boris_step(dTau/2, ...)`, which only reaches `t_final`, so the average of
p_{n-1/2} and p(t_final) sat at `t_final - dTau/4`. That fixed time offset made
`velo` and `momt` first order. Item 3 also limited momentum to first order, so
fixing either one alone shows little gain (see the tables below).

**Fix.** Use `boris_step(dTau, ...)`, giving p_{n+1/2}, so the centred average
lands on `t_final`. The `accl` diagnostic is now the centred difference
`(v_{n+1/2} - v_{n-1/2})/dTau`. The usage comment is updated.

**How it was checked** (scratch builds before and after, 2026-11-15; "before"
is the tree with fixes 1 to 3).

- *Against theory: momentum mid-pulse.* In a plane wave `p_perp + (q/m)*A_perp`
  is conserved, so the exact momentum at any time follows from the vector
  potential, with no dependence on the pusher. 20 electrons, gamma = 500,
  a0 = 10, `z0 = 2.5*Tau_FWHM`. The run is stopped at 45% of the window (where
  the field is on) at `--refine` 1, 2, 4, 8, so every grid ends at the same
  time. The table gives `max|p_perp - exact|/a0`, then the ratio between
  successive refinements (2 is first order, 4 is second order).

  | Case | `--halfstep` before | `--halfstep` after | Default path |
  |---|---|---|---|
  | Linear, 1 cycle | 2.2e-3 to 2.8e-4 (ratio 2.00) | 6.8e-6 to 1.1e-7 (ratio 4.00) | 4.5e-3 to 5.6e-4 (ratio 2.00) |
  | Linear, 4 cycles | 1.5e-3 to 1.9e-4 (2.00) | 4.3e-6 to 6.7e-8 (4.00) | 3.0e-3 to 3.7e-4 (2.00) |
  | Circular, 1 cycle | 1.6e-3 to 2.0e-4 (2.00) | 6.0e-6 to 9.3e-8 (4.00) | 3.2e-3 to 3.9e-4 (2.00) |
  | Circular, 4 cycles | 1.1e-3 to 1.3e-4 (2.00) | 3.6e-6 to 5.6e-8 (4.00) | 2.1e-3 to 2.6e-4 (2.00) |

  The error is second order only with the fix. The default path, which does not
  stagger, stays first order, as expected. This test also confirms the sign and
  normalisation of the Model 0 field from item 2 independently of the code.
- *Convergence (`convtest --order 2 --tol 0.25`, `--halfstep`).* Orders on
  posi, velo, momt:

  | Case | Before | After, refine 1/2/4 | After, refine 2/4/8 |
  |---|---|---|---|
  | Model 0, reference | 2.00, 1.00, 1.00 | 2.00, 2.00, 2.00 (PASS) | 2.00, 2.00, 2.00 (PASS) |
  | Model 2, `pola` 0 | 2.00, 1.00, 1.00 | 2.00, 2.00, 2.00 (PASS) | 2.00, 2.00, 2.00 (PASS) |
  | Model 2, `pola` 1 | 2.00, 1.00, 1.00 | 2.00, 2.00, 2.00 (PASS) | 2.00, 2.00, 2.00 (PASS) |

  `accl` and `egama` are ungated. After items 3 and 4 their errors are about
  1e-12 and 1e-9, which is round-off, so their observed orders (between -1.6 and
  1.4) mean nothing. I expected `accl` to reach order 2, and it did not, because
  of this floor.
- *No change to the default path.* Without `--halfstep`, old and new `bench_cpu`
  are `BIT-IDENTICAL` for Model 0 and for Model 2 with `pola` 0 and 1.
  Armadillo vs `bench_cpu` is still `BIT-IDENTICAL` for the Model 0 reference.
  `bench_gpu` was not touched.
- *Free drift.* `--a0 0 --halfstep`: momentum unchanged exactly, and position
  equals `x0 + v*T` to 1e-13 (positions of order 20).

**After pulling this change.** Only `--halfstep` output changes. Reference data
is not affected.

**Still first order.** LEADS, `bench_gpu`, and `bench_cpu` without `--halfstep`
start momentum at t0, not t0 - dTau/2, and report it half a step behind position.
They are first order wherever the field at the window ends is not negligible.
Staggering them as well would break bit-identity with existing reference data
unless `leads_dyna.cpp`, `bench_cpu` and `bench_gpu` all change together.
