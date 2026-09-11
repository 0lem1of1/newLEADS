/*
  Three-grid Richardson convergence-order tool for the ported Boris pusher.

  Kept deliberately separate from statediff.cpp: statediff's strict
  NP/nT-header-match guard exists to catch a real bug (comparing runs that
  aren't actually the same problem), and must never be weakened just to make
  refinement testing "easier" -- that's exactly what this tool is for
  instead. convtest only requires NP to match across the three inputs (same
  particles); nT/dTau are expected to differ, since these are meant to be
  the same physical run at coarse/medium/fine time resolution (bench_cpu
  --refine 1/N/M with matching --nsteps so every grid lands on the same
  final time).

  For three solutions at step sizes h, h/r, h/r^2 (r = refinement ratio):
      p = log( ||coarse-medium|| / ||medium-fine|| ) / log(r)
  This assumes a constant ratio r between successive grids -- cross-checked
  from the three files' own dTau headers rather than assumed to be 2.

  posi/velo/momt are the "gated" arrays: with --order X --tol Y, all three
  must satisfy |observed_order - X| <= Y for a pass (exit 0).

  egama/accl are reported but explicitly "ungated": both carry an O(dTau^2)
  Boris-invariant artifact that is expected behavior of a leapfrog-family
  integrator, not a defect -- gating on them would fail runs that are
  otherwise correct.

  Usage:
    convtest [--order X --tol Y] coarse.bin medium.bin fine.bin
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "../leads_state.h"

static double l2diff(const double* a, const double* b, long n)
{
    double sumsq = 0.0;
    for(long i = 0; i < n; i++) { double d = a[i]-b[i]; sumsq += d*d; }
    return sqrt(sumsq);
}

static void Usage(const char* argv0)
{
    fprintf(stderr,"Usage: %s [--order X --tol Y] coarse.bin medium.bin fine.bin\n",argv0);
}

int main(int argc, char** argv)
{
    double expect_order = NAN, order_tol = NAN;
    bool   gate = false;
    const char* paths[3] = {nullptr,nullptr,nullptr};
    int    npos = 0;

    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i],"--order") == 0 && i+1 < argc) { expect_order = atof(argv[++i]); }
        else if(strcmp(argv[i],"--tol") == 0 && i+1 < argc) { order_tol = atof(argv[++i]); }
        else if(npos < 3) { paths[npos++] = argv[i]; }
        else { Usage(argv[0]); return 1; }
    }
    if(npos != 3) { Usage(argv[0]); return 1; }
    if(!std::isnan(expect_order) && !std::isnan(order_tol)) gate = true;

    LeadsState c,m,f;
    if(!LeadsReadState(paths[0],c)) { fprintf(stderr,"Failed to read %s\n",paths[0]); return 1; }
    if(!LeadsReadState(paths[1],m)) { fprintf(stderr,"Failed to read %s\n",paths[1]); return 1; }
    if(!LeadsReadState(paths[2],f)) { fprintf(stderr,"Failed to read %s\n",paths[2]); return 1; }

    if(c.hd.NP != m.hd.NP || m.hd.NP != f.hd.NP)
    {
        fprintf(stderr,"NP mismatch across the three inputs (%ld, %ld, %ld) -- "
                        "these must be the same particle set\n",
                        (long)c.hd.NP,(long)m.hd.NP,(long)f.hd.NP);
        return 1;
    }
    const long NP = (long)c.hd.NP;

    // Refinement ratio cross-check from the files' own headers, not assumed.
    const double r_cm = c.hd.dTau / m.hd.dTau;
    const double r_mf = m.hd.dTau / f.hd.dTau;
    if(fabs(r_cm - r_mf) > 1e-6 * r_cm)
    {
        fprintf(stderr,"WARNING: refinement ratio coarse/medium (%.6g) != medium/fine (%.6g) -- "
                        "observed-order formula assumes a constant ratio; results below use "
                        "r=%.6g (coarse/medium)\n", r_cm, r_mf, r_cm);
    }
    const double r = r_cm;

    struct { const char* name; double* pc; double* pm; double* pf; long n; bool gated; } arrays[] = {
        {"posi",  c.posi,  m.posi,  f.posi,  3*NP, true},
        {"velo",  c.velo,  m.velo,  f.velo,  3*NP, true},
        {"momt",  c.momt,  m.momt,  f.momt,  3*NP, true},
        {"accl",  c.accl,  m.accl,  f.accl,  3*NP, false},
        {"egama", c.egama, m.egama, f.egama, NP,   false},
    };

    printf("Refinement ratio: coarse/medium=%.6g medium/fine=%.6g\n\n", r_cm, r_mf);
    printf("%-6s %14s %14s %10s %8s\n","array","||c-m||","||m-f||","order p","gated");

    bool all_pass = true;
    for(auto& arr : arrays)
    {
        double d_cm = l2diff(arr.pc, arr.pm, arr.n);
        double d_mf = l2diff(arr.pm, arr.pf, arr.n);
        double p = log(d_cm/d_mf) / log(r);

        printf("%-6s %14.6e %14.6e %10.4f %8s", arr.name, d_cm, d_mf, p, arr.gated ? "yes" : "no");
        if(!arr.gated)
            printf("  (ungated: known O(dTau^2) Boris-invariant artifact, not a bug)");
        printf("\n");

        if(gate && arr.gated)
        {
            bool pass = fabs(p - expect_order) <= order_tol;
            if(!pass) all_pass = false;
        }
    }

    if(gate)
    {
        printf("\n--order %.4g --tol %.4g on posi/velo/momt: %s\n",
               expect_order, order_tol, all_pass ? "PASS" : "FAIL");
    }

    LeadsFreeState(c); LeadsFreeState(m); LeadsFreeState(f);

    if(gate) return all_pass ? 0 : 1;
    return 0;
}
