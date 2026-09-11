/*
  Strict bit/tolerance comparator between two particle-state binary files
  (leads_state.h format) -- e.g. Armadillo's state_final_arma.bin against
  bench_cpu's or bench_gpu's state_final_scalar.bin / state_final_gpu.bin.

  Deliberately strict on header matching: NP/nT mismatches are a hard
  rejection, not a warning. Comparing two runs at different grid resolutions
  is what convtest.cpp is for -- do not weaken this guard to make that
  "easier", the two tools serve opposite goals.

  Usage:
    statediff [--tol X] fileA fileB

  Without --tol, this only reports (exit 0 always). With --tol X, exits 0 if
  the worst relative error across all five arrays is <= X, else 1.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "../leads_state.h"

// Floor on the denominator of a relative-error ratio, so a comparison isn't
// blown up by a value that's legitimately near zero (e.g. accl components
// that are ~0 by physics, not by bug).
static const double TINY = 1e-12;

struct DiffResult
{
    double max_abs;
    long   max_abs_idx;
    double max_rel;
    long   max_rel_idx;
    double l2;
};

static DiffResult compare(const double* a, const double* b, long n)
{
    DiffResult r{0,0,0,0,0};
    double sumsq = 0.0;
    for(long i = 0; i < n; i++)
    {
        double d = fabs(a[i]-b[i]);
        if(d > r.max_abs) { r.max_abs = d; r.max_abs_idx = i; }

        double denom = fabs(a[i]) > TINY ? fabs(a[i]) : TINY;
        double rel = d/denom;
        if(rel > r.max_rel) { r.max_rel = rel; r.max_rel_idx = i; }

        sumsq += d*d;
    }
    r.l2 = sqrt(sumsq);
    return r;
}

static void Usage(const char* argv0)
{
    fprintf(stderr,"Usage: %s [--tol X] fileA fileB\n",argv0);
}

int main(int argc, char** argv)
{
    double tol = -1.0;
    bool   have_tol = false;
    const char* pathA = nullptr;
    const char* pathB = nullptr;

    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i],"--tol") == 0 && i+1 < argc) { tol = atof(argv[++i]); have_tol = true; }
        else if(!pathA) pathA = argv[i];
        else if(!pathB) pathB = argv[i];
        else { Usage(argv[0]); return 1; }
    }
    if(!pathA || !pathB) { Usage(argv[0]); return 1; }

    LeadsState a,b;
    if(!LeadsReadState(pathA,a)) { fprintf(stderr,"Failed to read %s\n",pathA); return 1; }
    if(!LeadsReadState(pathB,b)) { fprintf(stderr,"Failed to read %s\n",pathB); return 1; }

    if(a.hd.NP != b.hd.NP || a.hd.nT != b.hd.nT)
    {
        fprintf(stderr,"HEADER MISMATCH: NP %ld vs %ld, nT %ld vs %ld -- refusing to compare "
                        "(use convtest for different grid resolutions)\n",
                        (long)a.hd.NP,(long)b.hd.NP,(long)a.hd.nT,(long)b.hd.nT);
        return 1;
    }

    const long NP = (long)a.hd.NP;

    struct { const char* name; double* pa; double* pb; long n; } arrays[] = {
        {"posi",  a.posi,  b.posi,  3*NP},
        {"velo",  a.velo,  b.velo,  3*NP},
        {"momt",  a.momt,  b.momt,  3*NP},
        {"accl",  a.accl,  b.accl,  3*NP},
        {"egama", a.egama, b.egama, NP},
    };

    printf("%-6s %14s %10s %14s %10s %14s\n","array","max_abs","@idx","max_rel","@idx","L2");
    double worst_rel = 0.0;
    for(auto& arr : arrays)
    {
        DiffResult r = compare(arr.pa, arr.pb, arr.n);
        printf("%-6s %14.6e %10ld %14.6e %10ld %14.6e\n",
               arr.name, r.max_abs, r.max_abs_idx, r.max_rel, r.max_rel_idx, r.l2);
        if(r.max_rel > worst_rel) worst_rel = r.max_rel;
    }

    if(worst_rel == 0.0) printf("BIT-IDENTICAL\n");
    else printf("Worst relative error across all arrays: %.6e\n", worst_rel);

    LeadsFreeState(a);
    LeadsFreeState(b);

    if(have_tol) return worst_rel <= tol ? 0 : 1;
    return 0;
}
