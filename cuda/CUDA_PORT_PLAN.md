# CUDA port plan — Model-0 Boris pusher

Companion to `KNOWN_ISSUES.md`. That file tracks bugs; this one tracks the
CPU -> CUDA port itself: what's done, what's decided, what's open.

Status key: `[ ]` not started · `[~]` in progress · `[x]` done

---

## 1. Scope

What's in / out of this port. (The CPU port already drew these lines —
confirm they still hold for CUDA or note where they change.)

- Model: **0 only** (Plain Wave), via `laser_profile<0>()` in `leads_boris.cuh`
- Radiation reaction: **off** (`RRFLAG == 0` branch only, i.e. `boris_step()`
  as already ported)
- Reference implementation: `leads_traj.cpp::CalForceMovePart` +
  `leads_laser.cpp::LaserProfile` case 0, via the scalar CPU port in
  `cuda/bench_cpu.cpp`
- [x] Anything else in scope / explicitly deferred?
  Nothing else identified. Trajectory history (full per-step cubes) is
  explicitly deferred — see section 3.

---

## 2. Target: what `bench_gpu.cu` needs to do

Match `bench_cpu.cpp`'s contract so `statediff`/`convtest` work unmodified
against it.

- Inputs: `state_init.bin` (`leads_state.h`), `params.bin` (`leads_params.h`)
- Output: same `LeadsState` binary format, written from the device result
- CLI surface to replicate (or explicitly drop) from `bench_cpu.cpp`:
  `--a0`, `--refine`, `--halfstep`, `--nsteps`
  - [x] `--halfstep` in the GPU port too, or defer until
        `KNOWN_ISSUES.md`'s convergence bug is resolved on CPU first?
    Decided: dropped from `bench_gpu` entirely — Armadillo has no halfstep
    resync either, so the GPU port only needs to match Armadillo, not
    replicate a CPU-only diagnostic. `bench_gpu`'s CLI is
    `--a0 X --refine N --nsteps N`, plus `--block N` (new, GPU-only).
- [x] Single-particle-per-thread, or something else (this pusher has no
      inter-particle coupling — Model 0's field is a function of `t,x,y,z`
      only, so this is embarrassingly parallel over particles)?
    Implemented: one thread per particle in `BorisKernel` (`cuda/bench_gpu.cu`).

---

## 3. Data layout on device

`leads_state.h` already stores `posi/velo/momt/accl` column-major
(`[x0..xN-1, y0..yN-1, z0..zN-1]`) — i.e. already SoA, which is what you want
for coalesced access per-component across particles.

- [x] Copy `LeadsState`'s flat arrays to device as-is (one `cudaMalloc` +
      `cudaMemcpy` per array), or repack first?
    Implemented: no repacking. `posi/velo/momt/accl` copied straight
    into `d_posi/d_velo/d_momt/d_accl` with the same `[0*NP+part,
    1*NP+part, 2*NP+part]` indexing `bench_cpu.cpp` uses.
- [x] `LeadsParams` (200-byte POD, no pointers) — `__constant__` memory
      (same params for every particle/thread) or plain global + pass by
      value?
    Went with `__constant__ LeadsParams d_params` + one `cudaMemcpyToSymbol`
    before launch, rather than passing it as a kernel argument. Note: on
    this toolchain (CUDA 13.2 / sm_86) a struct kernel argument is itself
    placed in constant memory by the compiler, so "pass by value" and
    explicit `__constant__` end up nearly equivalent here — `__constant__`
    was chosen because it's the more portable/explicit choice across
    toolchains, not because it measurably beat passing by value. A raw
    pointer to a device-side copy (dereferenced per-thread from global/L2)
    was **not** used since it would only add a layer of indirection for no
    benefit given the struct's size.
- [x] Per-particle history (`TR`, full trajectory cubes) — in scope for the
      GPU port, or is this harness final-state-only like `bench_cpu`
      currently is? (`leads_dyna.cpp`'s `StoreTrajectories` flag exists on
      the Armadillo side for this.)
    Decided: final-state-only, matching `bench_cpu`. `TR` is still built
    and written into the output file host-side (statediff/convtest need
    it) but is never copied to the device — the kernel recomputes
    `t = -t_shift + i*dTau` inline instead of indexing a device array,
    avoiding an `nT`-sized transfer for a value every thread can derive
    for free.

---

## 4. Kernel design

- [x] Grid/block sizing strategy for `NP` particles:
    Implemented: 256 threads/block default (multiple of the 32-thread warp
    size), grid size `ceil(NP/block)`, with an in-kernel `if(part >= NP)
    return` guard for the remainder. Block size is overridable via
    `--block N` for benchmarking.

    **Found a real limit here, not just a style choice:** `--block 1024`
    fails at launch with `too many resources requested for launch`.
    `nvcc -Xptxas -v` reports `BorisKernel` uses **68 registers/thread**;
    68 * 1024 = 69,632 exceeds `sm_86`'s 65,536-register file per SM, so
    1024 threads/block can never be resident regardless of occupancy
    settings. 32/128/256/512 were all tested and produced bit-identical
    output to each other (see section 6) — the register ceiling is a hard
    launch-configuration limit, not a correctness bug, but it means
    block sizes near 1024 need `--maxrregcount` or `launch_bounds` if ever
    required (not needed at 256).
- [x] Where does the timestep loop live — inside the kernel (one launch,
      `nsteps` iterations per thread), or one kernel launch per step (`nT`
      launches)? The former avoids launch overhead but forces the loop
      bound (`--refine`/`--nsteps`) to be a kernel-time argument, not a
      host-side loop like `bench_cpu.cpp` has.
    Implemented: loop lives inside `BorisKernel`, one launch total. Revisit
    only if per-step diagnostics (e.g. writing trajectory history back to
    the host) are needed later.
- [x] `boris_step()` and `laser_profile<0>()` reused verbatim from
      `leads_boris.cuh` — confirm no changes needed once compiled with
      `nvcc` (the header already guards this via `LEADS_HD`/`__CUDACC__`).
    Confirmed: `cuda/bench_gpu.cu` includes `leads_boris.cuh` unchanged and
    compiles clean under `nvcc -std=c++17 -O3 -gencode arch=compute_86,code=sm_86`.
- [x] Double precision throughout (matches `bench_cpu.cpp` and the
      Armadillo reference) — any reason to consider mixed precision, or is
      that explicitly out of scope given the machine-precision-match goal
      described in `leads_boris.cuh`'s header comment?
    Confirmed: FP64 throughout, no mixed precision. Consumer GPUs (this
    RTX 3050 included) run FP64 at a fraction of their FP32 throughput,
    but machine-precision matching against Armadillo is the actual goal
    here, not raw FLOPs — see section 7 for what that costs in practice.

---

## 5. Build

- [x] Add a `cuda/Makefile` (or extend the root `Makefile`) with a
      `bench_gpu` target compiling `bench_gpu.cu` via `nvcc`
    Added `cuda/Makefile` with `bench_cpu`/`bench_gpu`/`convtest`/`statediff`
    targets. `bench_gpu`: `nvcc -std=c++17 -O3 -gencode arch=$(ARCH)`.
- [x] Target GPU architecture(s) (`-arch=sm_XX` / `-gencode`) —
      what hardware is this expected to run on?
    Default `ARCH=compute_86,code=sm_86` — this dev machine's GPU
    (`nvidia-smi`: RTX 3050 Laptop, compute capability 8.6). Overridable
    per-build: `make bench_gpu ARCH=compute_75,code=sm_75`. No multi-arch
    fat binary yet — add more `-gencode` entries to `NVCCFLAGS` if this
    needs to run on other hardware later.

---

## 6. Validation plan

Reuse the tools already built for the CPU port — they diff two
`LeadsState` binaries, so they don't care whether the second one came from
`bench_cpu` or `bench_gpu`.

- [x] `statediff state_final_arma.bin bench_gpu_out.bin` — GPU vs Armadillo
      reference, same tolerance bar the CPU port had to clear
    Run on `TestCase/bench/` (NP=5000, nT=8001): worst relative error
    1.626e-01 (in `accl` only, at a single index where the finite-difference
    acceleration crosses zero — a known artifact, see `accl`'s row). posi/
    velo/momt/egama all at ~1e-6 relative / ~1e-9 absolute or better.
    **This is the same error, to the last printed digit, as `bench_cpu` vs
    Armadillo already produces** — the GPU port doesn't add any discrepancy
    of its own.
- [x] `statediff` CPU-output vs GPU-output directly — should match at
      (near-)bit level since both run the identical `leads_boris.cuh` math
      in the same precision
    Run on the same case: max abs error ~1.2e-11 (momt/egama), ~1.4e-14
    (posi), all near machine epsilon for accumulated FP64 error over 8001
    steps. `accl`'s 1.17e-3 relative error is the same zero-crossing
    artifact as above, not new. `--a0 0` free-drift case is **bit-identical**
    between CPU and GPU (checked directly, no tolerance needed).
- [x] `convtest` 3-grid Richardson test on the GPU path — expect the same
      order-1 result as CPU until `KNOWN_ISSUES.md`'s halfstep bug is
      fixed (this is a pusher/init-condition issue, not a CPU-vs-GPU one,
      so the GPU port should reproduce it, not diverge from it)
    Ran it anyway (cheap, and it's the strongest evidence the port is
    correct, not just "close enough"): observed order ~1.00 for
    posi/velo/momt, i.e. **exactly the CPU behavior**, including which
    rows are ungated (`accl`/`egama`, same known O(dTau^2) artifact). No
    GPU-specific divergence. This is not a bug to fix here — it's the same
    open item tracked in `KNOWN_ISSUES.md`, and fixing it there should
    make `bench_gpu` converge at order 2 automatically since both drivers
    share `leads_boris.cuh`.
- [x] New failure modes to check for that don't exist on CPU: race
      conditions (shouldn't apply — no cross-particle writes), NaN/inf from
      uninitialized device memory, block-size-dependent results (would
      indicate a bug — this problem has no shared/reduction state)
    Tested `--block 32/128/256/512`: all bit-identical to each other (no
    block-size dependence, as expected with zero shared state). `--block
    1024` **fails to launch** — not a correctness bug, a register-pressure
    limit, see section 4's grid/block finding. No NaN/inf observed in any
    run.

---

## 7. Performance

`cuda/bench_cpu` is the baseline. What's the actual goal here — is this port
about correctness/feasibility first, or is there a throughput target?

- [x] Baseline: `bench_cpu` wall-clock for `NP=5000`, `nT=8001`
    `PusherLoopTime_sec := 3.197` (`1.25e7` particle-steps/sec)
- [x] GPU wall-clock, same problem size
    Kernel only (`cudaEvent` timed): `PusherLoopTime_sec := 0.258`
    (`1.55e8` particle-steps/sec) — **~12.4x** faster than CPU.
    End-to-end including `cudaMalloc`/H2D/D2H/event overhead:
    `WallTime_sec := 0.449` — **~7.1x** faster than CPU.
- [x] Bottleneck expected to be launch overhead / memory transfer /
      compute, given how cheap `boris_step()` is per particle per step?
    At this problem size the gap between kernel time (0.258s) and wall
    time (0.449s) is transfer/setup overhead (~0.19s), not compute — and
    NP=5000 only fills 20 blocks of 256 threads, nowhere near enough to
    saturate this GPU's SMs. Both effects point the same way: this
    problem size is too small to show the GPU's real ceiling. Re-measure
    at larger NP (more particles, not more steps) before treating 12.4x/
    7.1x as representative — expect the kernel-only number to hold or
    improve with more particles (better occupancy), and the wall-time
    number to close in on it as transfer cost amortizes.

---

## 8. Open questions / decisions log

Running log of anything decided along the way that isn't obvious from the
code — mirrors how `KNOWN_ISSUES.md` documents the halfstep investigation.
Add dated entries as you go:

- **2026-09-13** — Initial `bench_gpu.cu` written and validated end-to-end
  against `TestCase/bench/` (NP=5000, nT=8001, Model 0). All checks in
  section 6 pass at parity with the already-validated CPU port. `cuda/Makefile`
  added. No open GPU-specific bugs; the only open issue is the pre-existing
  order-1 convergence bug tracked in `KNOWN_ISSUES.md`, which both drivers
  now share identically since they run the same `leads_boris.cuh` code.
- **2026-09-13** — Register pressure caps usable block size on `sm_86`:
  68 regs/thread * 1024 threads = 69,632 > the 65,536-register file per SM,
  so `--block 1024` fails to launch outright (963 is the actual per-block
  ceiling at 68 regs/thread, but block sizes are normally kept a power of
  two / multiple of 32, so 512 is the largest sane value below it). Not a
  problem at the default 256, noted here so it isn't rediscovered by
  surprise on a GPU with a smaller register file.

---

## 9. Milestones

- [x] `bench_gpu.cu` compiles and runs on a trivial case (`--a0 0` free
      drift, matches `bench_cpu --a0 0` via `statediff`) — bit-identical
- [x] Matches `bench_cpu` output on a real Model-0 case within tolerance
- [x] Matches Armadillo `state_final_arma.bin` within the same tolerance
      the CPU port achieved
- [x] `convtest` convergence order on GPU matches CPU (including
      reproducing the known halfstep issue, not fixing it here)
- [x] Performance numbers recorded (section 7) — re-measure at larger NP
      before treating the current 12.4x/7.1x figures as final
- [x] `cuda/Makefile` target added, build documented

Ported and validated. What's left is tuning (section 7's re-measurement at
larger NP) and the shared numerics bug in `KNOWN_ISSUES.md` — not new GPU
work.
