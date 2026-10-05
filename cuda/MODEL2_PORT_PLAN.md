# Model 2 port plan — design decisions

Companion to `CUDA_PORT_PLAN.md` (Model 0's port doc — that one's done and
validated). Sections 2-5 below are decided (2026-10-05); the harness and
port follow.

Status key: `[ ]` not started · `[~]` in progress · `[x]` decided/done

---

## 1. What Model 2 actually is

From `leads_laser.cpp` case 2: a tightly focused 3D pulse
(complex-source-point construction), not a plane wave like Model 0. Field
comes out of `real()`/`imag()` of ~12 `std::complex<double>` intermediates
(`Rc, tc, tau, p0, phase0, P_0, Zc, f, g, hf`, then the
`dummyEx/Ey/Ez/Bx/By/Bz` outputs), using complex `sqrt`, `exp`, complex
division, and one branch selection: `Rc = imag(Rc)<0 ? -Rc : Rc`. No
`E_Const`/`M_Const` normalization (same as Model 0).

Facts verified against the source:
- `Rc = sqrt(x² + y² + (z + i·z0)²)`, `z0 = k*zr` (a complex-source-point
  distance). x, y, z are the normalized `x1,y1,z1` (no `L_Const`).
- `pola` is `int` on the Armadillo side too (`leads.h:53`,
  `pola = int(in.Polarization)`), so non-integer ellipticity is not
  expressible even in the reference.
- Input defaults: `Wavelength=800nm`, `BeamWaist=3µm` (~3.75λ),
  `Tau_FWHM=4*2π`, `a0=10`, `PulseModel=0`, `Polarization=0`.
  `zr = k*w0²/2`, so `z0 = k*zr ≈ 275`.
- `a0` enters only via `E0 = B0 = a0` → `p0`, a linear prefactor on every
  field component, so `a0=0` gives exactly 0 field provided no NaN/Inf
  arises elsewhere (to be confirmed in the harness).

---

## 2. Complex-number strategy

- [x] **Q1 Complex type.** Dual path: `std::complex` under g++,
      `thrust::complex` under nvcc, split with `#ifdef __CUDACC__` like
      `LEADS_HD`. **Before** writing `laser_profile<2>()`, verify with a
      standalone harness (`cuda/cplx_harness.cu`) that evaluates the Model 2
      intermediates on host and device at sampled points, including the
      `z=0` plane. The `Rc` branch is the main risk. Harness covers complex
      division as well as `sqrt`/`exp`.
- [x] **Q2 Agreement bar.** Not bit-for-bit. Two tiers:
    - *Field level (harness):* ulp-level differences are fine; gross
      difference or sign flip blocks the port. The 1e-13 threshold floated
      earlier was a guess — set it after seeing the real error
      distribution.
    - *End to end:* `statediff`'s existing tolerance, unchanged. A
      CPU-vs-GPU mismatch that passes the harness is documented as
      accumulated rounding (like `KNOWN_ISSUES.md`), not a blocker.
- [x] **Q3 Fallback** (only if the harness fails), in order:
    1. Rerun with FMA contraction off (`-fmad=false`, `-ffp-contract=off`)
       to separate library effects from contraction.
    2. If a real mismatch remains, hand-write only the offending operation
       as `__host__ __device__`.
    3. If still unresolved, use the hand-written version on both CPU and GPU.
    4. A looser tolerance is never acceptable for field-level mismatches.

---

## 3. The `pola` field

- [x] **Q4 Type.** Keep `pola` as `int32_t`; Model 2 uses `(double)p.pola`.
      Covers linear (0) and circular (±1). Verified: Armadillo side is also
      `int pola`, so this matches the reference exactly.
- [x] **Q5 Non-integer ellipticity.** Not needed now. If added later,
      overlay a `double ellipticity` on one `reserved[]` slot with a union
      so existing `params.bin` files stay valid.

---

## 4. Model dispatch in `bench_cpu`/`bench_gpu`

- [x] **Q6 Dispatch.** Template driver and kernel on the model
      (`run_cpu<int Model>`, `kernel<int Model>`) with a single
      `switch(p.model)` before the loop/launch. Cases 0 and 2 only; any
      other value errors out. No per-step branching, no function pointers.
      Model 0's instantiation must stay unchanged — confirm by rerunning
      its existing `statediff`.
- [x] **Q7 CLI.** Model-agnostic. Model chosen only by `p.model` from
      `params.bin`, no `--model` flag. To check: `--halfstep` lives in the
      integrator (not in `laser_profile<0>`), and `--a0 0` gives exactly
      zero field for Model 2 (no NaN/Inf).

---

## 5. Reference test case

- [x] **Q8 Physical parameters.** Reuse Model 0's `Tau_FWHM` and `a0`.
      Model 0's `BeamWaist` is already 3 µm (~3.75λ), which is in the
      "few wavelengths" range, so keep it and record it. `a0 = 0` only for
      the free-drift check.
- [x] **Q9 Ellipticity values.** `pola = 0` first as a linear sanity check
      (not a numerical match to Model 0, since Model 2 is focused), then
      `pola = 1` as the main circular reference, optionally `pola = -1` for
      handedness. Separate Armadillo reference per value; field-level
      harness at all of them.
- [x] **Q10 Particles.**
    - *Primary:* new `state_init.bin`, Model 0's particle count, placed in
      the focal region (~±`w0` transverse, ~±Rayleigh length axial,
      including particles at/near `z=0`), checked to overlap the pulse in
      time.
    - *Secondary:* Model 0's existing `state_init.bin` as far-field /
      finiteness check and for timing comparison.
    - `state_init.bin` is written in `leads_dyna.cpp:53`
      (`./bench/state_init.bin`), from whatever the Armadillo run
      initializes — so focal-region placement means changing the
      initial-distribution inputs, not the dump hook.

---

## 6. Open questions / decisions log

- **2026-10-05:** §§2-5 decided. Dual-path complex type gated by a
  host-vs-device harness (two-tier agreement bar, staged fallback). `pola`
  stays `int32_t`. Drivers templated on model with one `switch`. CLI stays
  model-agnostic. Reference runs: `pola` 0, 1 (optionally -1), `w0` of a
  few wavelengths, focal-region particles plus Model 0's file as a
  secondary check.

- **2026-10-05 — harness result (`cplx_harness`, 1e6 points, pola 0/1/-1,
  focal-region + z=0 plane + on-axis + origin subsets, `w0`=23.56, `z0`=277.6
  in normalized units):** `std::complex` (host), `thrust::complex` (host) and
  `thrust::complex` (device) all agree to
  max|abs err|/global max|field| <= 8.3e-14, and each is within 6.2e-14 of a
  `std::complex<long double>` reference. No `Rc`-branch differences, no
  non-finite values, z=0 plane clean (<= 1.2e-14). `-fmad=false` /
  `-ffp-contract=off` changes nothing material (fallback step 1 not needed).
  `a0=0` gives exactly zero field everywhere, no NaN/Inf. **Caveat:** the
  *point-local* relative error has outliers (~1e3 of 7e5 points >1e-6, up to
  O(1)) in the far tails (|z|~200, field ~1e-13 of peak) where the formula
  cancels catastrophically; there the double-precision paths are equally far
  from the long-double truth, so this is the formula's conditioning, not a
  port defect. Pass/fail metric is therefore the global-scale one; tentative
  field-level gate for the port: <= 1e-12 of global max|field|.

- **2026-10-05 — first Model 2 smoke test** (far-field: Model 0's
  `state_init.bin`, params patched to model=2, pola 0 and 1, 5000 particles x
  8001 steps): CPU-vs-GPU `statediff` L2 ~1e-10 on momt/egama, ~3e-13 on
  posi (same size as Model 0's CPU-vs-GPU gap); `--a0 0` free drift is
  BIT-IDENTICAL CPU vs GPU; GPU 0.44 s vs CPU 10.6 s. (An earlier note
  here called these particles "far-field"/weak evidence -- wrong: the
  default bunch sits at z0=20, radius 1 in 1/k units, i.e. well inside the
  focal region (w0=23.6, z0=277.6), and the field moves momt by up to 2e-2
  vs a0=0 (L2 0.71), far above the 2e-10 CPU/GPU gap. So Model 0's
  `state_init.bin` already is a valid focal-region test and the planned
  separate focal-region particle file is not needed.)
- **2026-10-05 -- Armadillo reference + validation (Model 2).** References
  generated from a scratch copy of the repo with `PulseModel=2`,
  `Polarization`=0 and 1 (`./leads 4`, ~25 s each), stored in
  `TestCase/bench_model2_pola{0,1}/` (untracked, like the rest of
  `TestCase/`). Model 0's `Tau_FWHM`, `a0`, `BeamWaist` and particle setup
  unchanged. Results (NP=5000, 8001 steps):
    - CPU vs Armadillo, pola=0: posi L2 2.1e-7 (max abs 2.8e-9), momt L2
      5.4e-5 (max abs 1.3e-6), egama L2 6.9e-8. pola=1: 2.1e-7 / 4.6e-5 /
      8.4e-8. Model 0's own gap vs Armadillo is posi 2.2e-7 / momt 5.8e-5 /
      egama 1.1e-8, so Model 2 is at the same bar (~1e-9 abs posi, ~1e-6
      abs momt).
    - GPU vs Armadillo: identical to the CPU numbers to the printed digits;
      GPU-vs-CPU gap ~1e-10 L2 (smoke test above).
    - `convtest` (refine 1/2/4, no `--halfstep`, pola=1): posi 0.9994, velo
      1.0000, momt 1.0000 -- first order, same as Model 0 without
      `--halfstep` (see `KNOWN_ISSUES.md`). `--halfstep` not retested for
      Model 2.
    - Block-size sweep 32/64/128/256/512: BIT-IDENTICAL outputs. 1024 fails
      to launch ("too many resources requested"): `BorisKernel<2>` uses 102
      registers/thread (vs 68 for Model 0), so max block size for Model 2 is
      512 (512 is also ~2x slower than 256: 0.82 s vs 0.42 s).

Resolved assumptions:
- `statediff` has no built-in tolerance: it reports, and only gates when
  given `--tol X` (relative). The working bar is the Model 0 baseline
  quoted above.
- `--halfstep` lives in the `bench_cpu` driver loop, not in the laser, so it
  is model-agnostic (not retested for Model 2).
- Initial distribution comes from `leads_init.cpp` (GSL Gaussian radius
  `BunchRadius`, flat length `BunchLength`, centered at `z0`); default
  bunch is already inside Model 2's focal region.
- `pola=-1` (optional handedness check) was not run against Armadillo; the
  field harness shows it identical to `pola=1`.

---

## 7. Milestones

- [x] Sections 1-5 above decided
- [x] `thrust::complex` vs `std::complex` agreement confirmed (no fallback
      needed) via `cuda/cplx_harness.cu`, before writing `laser_profile<2>()`
- [x] Model-2 reference data generated (`state_init.bin`/`params.bin`/
      `state_final_arma.bin`), `TestCase/bench_model2_pola{0,1}/`
- [x] `laser_profile<2>()` compiles under both g++ and nvcc
      (`LeadsComplex` in `leads_boris.cuh`)
- [x] `bench_cpu`/`bench_gpu` dispatch to Model 2 correctly (generic lambda
      `run_cpu` in `bench_cpu.cpp`, `BorisKernel<Model>` + one `switch` in
      `bench_gpu.cu`; other models error out). Model 0 CPU output is
      bit-identical to its pre-refactor output.
- [x] Validated: free-drift (`--a0 0`), `statediff` vs Armadillo,
      `statediff` CPU-vs-GPU, `convtest` convergence order, block-size sweep
