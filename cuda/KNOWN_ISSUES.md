# Known issues

## `bench_cpu --halfstep`: velo/momt/egama stuck at 1st-order convergence

**Status:** unresolved, deferred. Not blocking the CPU/GPU port itself --
`--halfstep` is a diagnostic extra on top of the raw pusher, not required for
matching Armadillo.

**Symptom:** with `--halfstep` off, `convtest`'s 3-grid Richardson test shows
`posi`, `velo`, and `momt` all converging at observed order ~1.0 (halving the
timestep only halves the error, not quarters it) -- expected, since
`leads_init.cpp` sets position and momentum at the same instant instead of
staggering them the half-step apart a leapfrog scheme wants.

`--halfstep` is meant to fix this: stagger momentum backward by a half-step
before the main loop, run the loop as normal, then resynchronize momentum
back onto position's time at the end via a centered average of the two
half-steps bracketing it.

**What's been tried:**

1. First cut had the initial half-kick's *direction* backwards (pushed
   momentum forward in time instead of backward). Fixing the sign alone took
   `posi`'s observed order from 1.00 -> **1.89**, matching prior recollection
   of this fix almost exactly. Confirms the backward initial half-kick is
   correct.
2. `velo`/`momt`/`egama` did NOT follow -- stayed at order ~1.0 even after
   reworking the end-of-loop resync to do a proper extra forward half-kick
   (evaluating the field at the final position/time) followed by a true
   centered average `(p_{n-1/2}+p_{n+1/2})/2`, rather than the original,
   sloppier attempt (averaging two already-computed full-step values at the
   wrong loop index). The error *magnitude* dropped meaningfully from that
   rework (momt's coarse-vs-medium L2 diff: 3.98e-3 -> 7.78e-4), so something
   real improved, but the *convergence rate* did not.

**Why this is puzzling:** position is only touched by the initial half-kick
(the resync logic runs after position is already finalized for the step), so
position's jump to order 1.89 shows the initial backward half-kick is
genuinely correct. Momentum's own per-step update uses that same,
already-validated Boris kick math, so by the same leapfrog-staggering theory
that explains position's improvement, momentum's own natural (staggered)
values should also be 2nd-order accurate at their half-integer grid points --
and averaging two 2nd-order-accurate points bracketing the target time should
still be 2nd-order. On paper momentum should have improved alongside
position. It didn't.

**Open hypotheses, untested:**
- The end-of-loop resync has a subtler bug not yet found.
- Reverse-integrating the Boris rotation via a negative `dTau` argument to
  `boris_step()` is a reasonable approximation for position (confirmed above)
  but not exact enough for momentum specifically.
- The original (lost) implementation used a genuinely different resync
  technique than what's been reconstructed here from memory.

**Where the code stands right now:** `cuda/bench_cpu.cpp`'s `--halfstep`
path implements attempt #2 above (the "proper" resync). It's left in this
state -- not reverted -- since it's a strict improvement in error magnitude
over attempt #1 even though the order-2 goal for velo/momt isn't met yet.

**To pick this back up:** reproduce with
```
bench_cpu bench/state_init.bin bench/params.bin out.bin --refine N --nsteps <(nT_orig-1)*N> --halfstep
```
at N=1,2,4 (using nsteps so each run lands exactly on +t_shift) and feed the
three outputs to `convtest coarse medium fine --order 2 --tol 0.25`.
