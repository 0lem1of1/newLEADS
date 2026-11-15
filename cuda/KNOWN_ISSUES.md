# Known issues

## `bench_cpu --halfstep`: velo/momt were stuck at 1st-order convergence (RESOLVED)

**Status:** resolved. See `../BUGFIXES.md`, items 3 and 4.

**Symptom (old).** With `--halfstep` off, `convtest`'s 3-grid Richardson test
shows `posi`, `velo`, and `momt` all converging at observed order ~1.0, as
expected: `leads_init.cpp` sets position and momentum at the same instant
instead of staggering them half a step apart as a leapfrog scheme wants.
`--halfstep` staggers momentum backward before the main loop and resynchronizes
it at the end. That took `posi` to order 2, but `velo` and `momt` stayed at ~1.

**Cause.** Two independent errors, and both had to be fixed:

1. *The Boris rotation used gamma from the start of the step* instead of from
   `UM` (the momentum after the first half E-kick). That left an O(dTau)
   error in every rotation. Fixed in `leads_traj.cpp` and `leads_boris.cuh`
   (`BUGFIXES.md` item 3).
2. *The end-of-loop resync kicked by `dTau/2`, not `dTau`.* After the loop,
   momentum is p_{n-1/2}. A kick of `dTau/2` reaches `t_final`, not
   `t_final + dTau/2`, so averaging with p_{n-1/2} landed at `t_final - dTau/4`,
   a fixed offset that is first order. The kick is now `dTau`, giving p_{n+1/2},
   and the average lands on `t_final` (`BUGFIXES.md` item 4).

Either fix alone leaves momentum near order 1. `convtest` orders with
`--halfstep` on the Model 0 reference case (refine 1/2/4):

| Fixes applied | posi | velo | momt |
|---|---|---|---|
| Neither (as first recorded here) | 1.89 | about 1.0 | about 1.0 |
| Item 3 only | 2.00 | 0.999 | 0.999 |
| Item 4 only | 2.03 | 1.45 | 1.06 |
| Both | 2.00 | 2.00 | 2.00 |

The Model 2 results are in `BUGFIXES.md`.

**The earlier hypotheses.**
- *A subtler bug in the resync:* yes, the `dTau/2` kick above.
- *Reverse-integrating the rotation with a negative `dTau` is not exact enough:*
  no. Boris is time-reversible, and the backward initial kick was correct (it
  took `posi` to order 2.00).
- *The lost original used a different technique:* not needed.

**What is still first order.** Only the `--halfstep` path is staggered.
LEADS itself, `bench_gpu`, and `bench_cpu` without `--halfstep` still start
momentum at t0 rather than t0 - dTau/2, and report it half a step behind
position. Their results are therefore first order wherever the field at the ends
of the window is not negligible (for example the reference case with
`z0 = 20`). Changing that means staggering in `leads_dyna.cpp` and `bench_gpu`
too, and it would break bit-identity with the existing reference data unless
all three change together.

**Reproduce.**
```
bench_cpu bench/state_init.bin bench/params.bin out.bin --refine N --nsteps $((8000*N)) --halfstep
```
at N=1,2,4, then `convtest coarse medium fine --order 2 --tol 0.25`.
