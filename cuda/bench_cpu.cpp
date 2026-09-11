/*
  Scalar CPU driver for the ported Boris pusher (Model-0 laser only).

  Loads a particle-state snapshot (state_init.bin, from leads_state.h) and the
  matching laser/particle params sidecar (params.bin, from leads_params.h),
  runs the pusher on the CPU using leads_boris.cuh, and writes the resulting
  final state in the same binary format so it can be diffed against
  Armadillo's own state_final_arma.bin with statediff (or against
  three grids' worth of output with convtest).

  No Armadillo, no GSL — this is meant to be the first of two ports (CPU here,
  CUDA in bench_gpu.cu) that both reuse leads_boris.cuh's exact math.

  Usage:
    bench_cpu <state_init.bin> <params.bin> <output.bin>
              [--a0 X] [--refine N] [--halfstep] [--nsteps N]

    --a0 X       Override params.a0 (X=0 gives a free-drift sanity check).
    --refine N   Rebuild the time grid at N x the resolution, spanning the
                 same [-t_shift, t_shift] window read from state_init.bin.
                 Feeds convtest's 3-grid Richardson convergence test.
    --halfstep   Best-effort reconstruction of a previously-validated
                 diagnostic: stagger momentum half a step behind position at
                 the start (single momentum-only half-kick before the main
                 loop), then resynchronize at the end via the centered
                 average (p_{n-1/2}+p_{n+1/2})/2, rather than a naive
                 backward re-push (which was found to double the error).
                 Flagged here as best-effort since it was rebuilt from
                 memory -- sanity-check its effect via --refine/convtest
                 before trusting it quantitatively.
    --nsteps N   Run exactly N pusher steps instead of the state file's own
                 nT. Without this flag the default is nT (matching
                 leads_dyna.cpp's own `for(i=0;i<nT;i++)` loop bit-for-bit,
                 including its one-step endpoint overshoot -- needed for
                 exact comparison against state_final_arma.bin). With
                 --refine, pass --nsteps (nT_orig-1)*N to land exactly on
                 +t_shift instead of overshooting.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>
#include "../leads_state.h"
#include "../leads_params.h"
#include "leads_boris.cuh"

static void Usage(const char* argv0)
{
    fprintf(stderr,
        "Usage: %s <state_init.bin> <params.bin> <output.bin>\n"
        "          [--a0 X] [--refine N] [--halfstep] [--nsteps N]\n", argv0);
}

// gamma - pz: exact analytic invariant of motion for a charged particle in a
// plane EM wave (field depends only on eta = t - z), independent of a0,
// polarization, or envelope shape. Model 0 only.
static inline double inv_gmpz(double px, double py, double pz)
{
    double gamma = sqrt(1 + (((px*px) + (py*py)) + (pz*pz)));
    return gamma - pz;
}

int main(int argc, char** argv)
{
    if(argc < 4) { Usage(argv[0]); return 1; }

    const char* init_path   = argv[1];
    const char* params_path = argv[2];
    const char* out_path    = argv[3];

    double a0_override = NAN;
    int    refine       = 1;
    bool   halfstep      = false;
    long   nsteps_override = -1;

    for(int i = 4; i < argc; i++)
    {
        if(strcmp(argv[i],"--a0") == 0 && i+1 < argc)      { a0_override = atof(argv[++i]); }
        else if(strcmp(argv[i],"--refine") == 0 && i+1 < argc) { refine = atoi(argv[++i]); }
        else if(strcmp(argv[i],"--halfstep") == 0)         { halfstep = true; }
        else if(strcmp(argv[i],"--nsteps") == 0 && i+1 < argc) { nsteps_override = atol(argv[++i]); }
        else { fprintf(stderr,"Unrecognized argument: %s\n",argv[i]); Usage(argv[0]); return 1; }
    }
    if(refine < 1) { fprintf(stderr,"--refine must be >= 1\n"); return 1; }

    LeadsState st;
    if(!LeadsReadState(init_path,st)) { fprintf(stderr,"Failed to read %s\n",init_path); return 1; }

    LeadsParams p;
    if(!LeadsReadParams(params_path,p)) { fprintf(stderr,"Failed to read %s\n",params_path); return 1; }

    if(!std::isnan(a0_override)) p.a0 = a0_override;

    const long NP       = (long)st.hd.NP;
    const long nT_orig   = (long)st.hd.nT;
    const double dTau_orig = st.hd.dTau;
    const double t_shift    = st.hd.t_shift;

    // Rebuild the time grid at N x resolution over the same span, matching
    // linspace's convention (nT points span nT-1 equal intervals).
    long   nT   = (nT_orig - 1) * refine + 1;
    double dTau = dTau_orig / refine;

    std::vector<double> TR(nT);
    for(long i = 0; i < nT; i++) TR[i] = -t_shift + i*dTau;

    // Default: replicate leads_dyna.cpp's own loop bit-for-bit, including its
    // one-step endpoint overshoot (`for(i=0;i<nT;i++)`, nT steps over nT-1
    // intervals) -- required for an exact match against state_final_arma.bin.
    long steps_to_run = (nsteps_override >= 0) ? nsteps_override : nT;
    if(steps_to_run > nT) { fprintf(stderr,"--nsteps exceeds refined grid size (%ld)\n",nT); return 1; }

    std::vector<double> posi(st.posi, st.posi+3*NP);
    std::vector<double> velo(st.velo, st.velo+3*NP);
    std::vector<double> momt(st.momt, st.momt+3*NP);
    std::vector<double> accl(st.accl, st.accl+3*NP);
    std::vector<double> egama(st.egama, st.egama+NP);

    double max_dgmpz = 0.0;

    auto PusherStart = std::chrono::high_resolution_clock::now();

    for(long part = 0; part < NP; part++)
    {
        double x = posi[0*NP+part], y = posi[1*NP+part], z = posi[2*NP+part];
        double px = momt[0*NP+part], py = momt[1*NP+part], pz = momt[2*NP+part];
        double vx = velo[0*NP+part], vy = velo[1*NP+part], vz = velo[2*NP+part];
        double ax = accl[0*NP+part], ay = accl[1*NP+part], az = accl[2*NP+part];

        const double gmpz0 = inv_gmpz(px,py,pz);

        if(halfstep)
        {
            // Stagger momentum HALF A STEP BEHIND position: step it backward
            // to t=TR[0]-dTau/2 with a negative half-timestep, not forward.
            double Ex,Ey,Ez,Bx,By,Bz;
            laser_profile<0>(TR[0], x,y,z, p, Ex,Ey,Ez,Bx,By,Bz);
            double dummy_x=0, dummy_y=0, dummy_z=0, dummy_ax,dummy_ay,dummy_az;
            boris_step(-dTau/2.0, p.q_part, p.m_part,
                       dummy_x,dummy_y,dummy_z, px,py,pz, vx,vy,vz, dummy_ax,dummy_ay,dummy_az,
                       Ex,Ey,Ez,Bx,By,Bz);
        }

        double px_prev = px, py_prev = py, pz_prev = pz;
        double vx_prev = vx, vy_prev = vy, vz_prev = vz;

        for(long i = 0; i < steps_to_run; i++)
        {
            if(halfstep && i == steps_to_run - 1)
            {
                px_prev = px; py_prev = py; pz_prev = pz;
                vx_prev = vx; vy_prev = vy; vz_prev = vz;
            }

            double Ex,Ey,Ez,Bx,By,Bz;
            laser_profile<0>(TR[i], x,y,z, p, Ex,Ey,Ez,Bx,By,Bz);
            boris_step(dTau, p.q_part, p.m_part,
                       x,y,z, px,py,pz, vx,vy,vz, ax,ay,az,
                       Ex,Ey,Ez,Bx,By,Bz);
        }

        double gamma_final;
        if(halfstep)
        {
            // Centered average (p_{n-1/2}+p_{n+1/2})/2 resynchronizes momentum
            // to position's time, rather than a naive backward re-push.
            px = 0.5*(px_prev + px);
            py = 0.5*(py_prev + py);
            pz = 0.5*(pz_prev + pz);
            gamma_final = sqrt(1 + (((px*px)+(py*py))+(pz*pz)));
            vx = px/gamma_final; vy = py/gamma_final; vz = pz/gamma_final;
            ax = (vx - vx_prev)/dTau;
            ay = (vy - vy_prev)/dTau;
            az = (vz - vz_prev)/dTau;
        }
        else
        {
            gamma_final = sqrt(1 + (((px*px)+(py*py))+(pz*pz)));
        }

        posi[0*NP+part]=x;  posi[1*NP+part]=y;  posi[2*NP+part]=z;
        velo[0*NP+part]=vx; velo[1*NP+part]=vy; velo[2*NP+part]=vz;
        momt[0*NP+part]=px; momt[1*NP+part]=py; momt[2*NP+part]=pz;
        accl[0*NP+part]=ax; accl[1*NP+part]=ay; accl[2*NP+part]=az;
        egama[part] = gamma_final;

        double dgmpz = fabs(inv_gmpz(px,py,pz) - gmpz0);
        if(dgmpz > max_dgmpz) max_dgmpz = dgmpz;
    }

    auto PusherEnd = std::chrono::high_resolution_clock::now();
    double PusherSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(PusherEnd-PusherStart).count();

    double checksum = 0.0;
    for(long i = 0; i < 3*NP; i++) checksum += fabs(posi[i]);
    for(long i = 0; i < 3*NP; i++) checksum += fabs(momt[i]);

    bool wrote = LeadsWriteState(out_path, NP, nT, dTau, t_shift,
                                  TR.data(), posi.data(), velo.data(),
                                  momt.data(), accl.data(), egama.data());
    if(!wrote) { fprintf(stderr,"Failed to write %s\n",out_path); return 1; }

    printf("PusherLoopTime_sec := %g\n", PusherSeconds);
    printf("FinalStateChecksum := %g\n", checksum);
    printf("ParticleStepsPerSec := %g\n", (double)NP * (double)steps_to_run / PusherSeconds);
    if(p.model == 0)
        printf("MaxAbsDeltaInvGammaMinusPz := %g\n", max_dgmpz);
    printf("Wrote %s (NP=%ld nT=%ld dTau=%g steps_run=%ld)\n", out_path, NP, nT, dTau, steps_to_run);

    return 0;
}
