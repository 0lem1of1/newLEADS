/*
  CUDA driver for the ported Boris pusher (Model-0 laser only).

  Same contract as bench_cpu.cpp: loads state_init.bin/params.bin (leads_state.h
  / leads_params.h), runs the pusher, writes a final-state binary in the same
  format so it can be diffed against state_final_arma.bin or bench_cpu's own
  output with statediff, or fed to convtest for the 3-grid Richardson test.

  Decisions recorded in CUDA_PORT_PLAN.md:
    - No --halfstep here: the Armadillo reference doesn't stagger momentum
      either, so this port only needs to match Armadillo, not add a diagnostic
      Armadillo never had.
    - One CUDA thread per particle (Model 0's field depends only on t,x,y,z,
      so particles never interact -- embarrassingly parallel).
    - LeadsState's arrays are already SoA (leads_state.h's column-major
      layout), copied to device unchanged.
    - LeadsParams (200-byte POD, same for every particle) lives in
      __constant__ memory.
    - The whole timestep loop runs inside one kernel launch per particle
      thread, not one launch per step -- avoids nT launches for what is a
      very cheap per-step computation.

  Usage:
    bench_gpu <state_init.bin> <params.bin> <output.bin>
              [--a0 X] [--refine N] [--nsteps N] [--block N]

    --a0 X       Override params.a0 (X=0 gives a free-drift sanity check).
    --refine N   Rebuild the time grid at N x the resolution, spanning the
                 same [-t_shift, t_shift] window read from state_init.bin.
                 Feeds convtest's 3-grid Richardson convergence test.
    --nsteps N   Run exactly N pusher steps instead of the state file's own
                 nT. Without this flag the default is nT (matching
                 leads_dyna.cpp's own `for(i=0;i<nT;i++)` loop bit-for-bit,
                 including its one-step endpoint overshoot -- needed for
                 exact comparison against state_final_arma.bin). With
                 --refine, pass --nsteps (nT_orig-1)*N to land exactly on
                 +t_shift instead of overshooting.
    --block N    CUDA block size (threads per block). Default 256.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>
#include <cuda_runtime.h>
#include "../leads_state.h"
#include "../leads_params.h"
#include "leads_boris.cuh"

#define CUDA_CHECK(call) do { \
    cudaError_t err__ = (call); \
    if(err__ != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err__)); \
        exit(1); \
    } \
} while(0)

__constant__ LeadsParams d_params;

static void Usage(const char* argv0)
{
    fprintf(stderr,
        "Usage: %s <state_init.bin> <params.bin> <output.bin>\n"
        "          [--a0 X] [--refine N] [--nsteps N] [--block N]\n", argv0);
}

// gamma - pz: exact analytic invariant of motion for a charged particle in a
// plane EM wave (field depends only on eta = t - z), independent of a0,
// polarization, or envelope shape. Model 0 only.
LEADS_HD static inline double inv_gmpz(double px, double py, double pz)
{
    double gamma = sqrt(1 + (((px*px) + (py*py)) + (pz*pz)));
    return gamma - pz;
}

// One thread per particle. Reads/writes the SoA arrays with stride NP,
// exactly like bench_cpu.cpp's indexing (posi[0*NP+part], posi[1*NP+part], ...).
__global__ void BorisKernel(long NP, long steps_to_run,
                             double t_shift, double dTau,
                             double* posi, double* velo,
                             double* momt, double* accl,
                             double* egama, double* dgmpz_out)
{
    long part = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if(part >= NP) return;

    double x = posi[0*NP+part], y = posi[1*NP+part], z = posi[2*NP+part];
    double px = momt[0*NP+part], py = momt[1*NP+part], pz = momt[2*NP+part];
    double vx = velo[0*NP+part], vy = velo[1*NP+part], vz = velo[2*NP+part];
    double ax = accl[0*NP+part], ay = accl[1*NP+part], az = accl[2*NP+part];

    const double gmpz0 = inv_gmpz(px,py,pz);

    for(long i = 0; i < steps_to_run; i++)
    {
        const double t = -t_shift + i*dTau;
        double Ex,Ey,Ez,Bx,By,Bz;
        laser_profile<0>(t, x,y,z, d_params, Ex,Ey,Ez,Bx,By,Bz);
        boris_step(dTau, d_params.q_part, d_params.m_part,
                   x,y,z, px,py,pz, vx,vy,vz, ax,ay,az,
                   Ex,Ey,Ez,Bx,By,Bz);
    }

    const double gamma_final = sqrt(1 + (((px*px)+(py*py))+(pz*pz)));

    posi[0*NP+part]=x;  posi[1*NP+part]=y;  posi[2*NP+part]=z;
    velo[0*NP+part]=vx; velo[1*NP+part]=vy; velo[2*NP+part]=vz;
    momt[0*NP+part]=px; momt[1*NP+part]=py; momt[2*NP+part]=pz;
    accl[0*NP+part]=ax; accl[1*NP+part]=ay; accl[2*NP+part]=az;
    egama[part] = gamma_final;

    dgmpz_out[part] = fabs(inv_gmpz(px,py,pz) - gmpz0);
}

int main(int argc, char** argv)
{
    if(argc < 4) { Usage(argv[0]); return 1; }

    const char* init_path   = argv[1];
    const char* params_path = argv[2];
    const char* out_path    = argv[3];

    double a0_override = NAN;
    int    refine       = 1;
    long   nsteps_override = -1;
    int    block_size    = 256;

    for(int i = 4; i < argc; i++)
    {
        if(strcmp(argv[i],"--a0") == 0 && i+1 < argc)      { a0_override = atof(argv[++i]); }
        else if(strcmp(argv[i],"--refine") == 0 && i+1 < argc) { refine = atoi(argv[++i]); }
        else if(strcmp(argv[i],"--nsteps") == 0 && i+1 < argc) { nsteps_override = atol(argv[++i]); }
        else if(strcmp(argv[i],"--block") == 0 && i+1 < argc)  { block_size = atoi(argv[++i]); }
        else { fprintf(stderr,"Unrecognized argument: %s\n",argv[i]); Usage(argv[0]); return 1; }
    }
    if(refine < 1) { fprintf(stderr,"--refine must be >= 1\n"); return 1; }
    if(block_size < 1) { fprintf(stderr,"--block must be >= 1\n"); return 1; }

    LeadsState st;
    if(!LeadsReadState(init_path,st)) { fprintf(stderr,"Failed to read %s\n",init_path); return 1; }

    LeadsParams p;
    if(!LeadsReadParams(params_path,p)) { fprintf(stderr,"Failed to read %s\n",params_path); return 1; }

    if(!std::isnan(a0_override)) p.a0 = a0_override;

    const long NP        = (long)st.hd.NP;
    const long nT_orig    = (long)st.hd.nT;
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
    std::vector<double> dgmpz(NP);

    auto WallStart = std::chrono::high_resolution_clock::now();

    CUDA_CHECK(cudaMemcpyToSymbol(d_params, &p, sizeof(LeadsParams)));

    double *d_posi, *d_velo, *d_momt, *d_accl, *d_egama, *d_dgmpz;
    const size_t n3_bytes = 3*(size_t)NP*sizeof(double);
    const size_t n1_bytes = (size_t)NP*sizeof(double);
    CUDA_CHECK(cudaMalloc(&d_posi, n3_bytes));
    CUDA_CHECK(cudaMalloc(&d_velo, n3_bytes));
    CUDA_CHECK(cudaMalloc(&d_momt, n3_bytes));
    CUDA_CHECK(cudaMalloc(&d_accl, n3_bytes));
    CUDA_CHECK(cudaMalloc(&d_egama, n1_bytes));
    CUDA_CHECK(cudaMalloc(&d_dgmpz, n1_bytes));

    CUDA_CHECK(cudaMemcpy(d_posi, posi.data(), n3_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_velo, velo.data(), n3_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_momt, momt.data(), n3_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_accl, accl.data(), n3_bytes, cudaMemcpyHostToDevice));

    cudaEvent_t KernelStart, KernelEnd;
    CUDA_CHECK(cudaEventCreate(&KernelStart));
    CUDA_CHECK(cudaEventCreate(&KernelEnd));

    const int grid_size = (int)((NP + block_size - 1) / block_size);

    CUDA_CHECK(cudaEventRecord(KernelStart));
    BorisKernel<<<grid_size, block_size>>>(NP, steps_to_run, t_shift, dTau,
                                            d_posi, d_velo, d_momt, d_accl,
                                            d_egama, d_dgmpz);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaEventRecord(KernelEnd));
    CUDA_CHECK(cudaEventSynchronize(KernelEnd));

    float KernelMs = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&KernelMs, KernelStart, KernelEnd));

    CUDA_CHECK(cudaMemcpy(posi.data(), d_posi, n3_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(velo.data(), d_velo, n3_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(momt.data(), d_momt, n3_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(accl.data(), d_accl, n3_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(egama.data(), d_egama, n1_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dgmpz.data(), d_dgmpz, n1_bytes, cudaMemcpyDeviceToHost));

    cudaFree(d_posi); cudaFree(d_velo); cudaFree(d_momt);
    cudaFree(d_accl); cudaFree(d_egama); cudaFree(d_dgmpz);
    cudaEventDestroy(KernelStart); cudaEventDestroy(KernelEnd);

    auto WallEnd = std::chrono::high_resolution_clock::now();
    double WallSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(WallEnd-WallStart).count();

    double max_dgmpz = 0.0;
    for(long part = 0; part < NP; part++) if(dgmpz[part] > max_dgmpz) max_dgmpz = dgmpz[part];

    double checksum = 0.0;
    for(long i = 0; i < 3*NP; i++) checksum += fabs(posi[i]);
    for(long i = 0; i < 3*NP; i++) checksum += fabs(momt[i]);

    bool wrote = LeadsWriteState(out_path, NP, nT, dTau, t_shift,
                                  TR.data(), posi.data(), velo.data(),
                                  momt.data(), accl.data(), egama.data());
    if(!wrote) { fprintf(stderr,"Failed to write %s\n",out_path); return 1; }

    const double KernelSeconds = KernelMs / 1000.0;
    printf("PusherLoopTime_sec := %g\n", KernelSeconds);
    printf("WallTime_sec := %g\n", WallSeconds);
    printf("FinalStateChecksum := %g\n", checksum);
    printf("ParticleStepsPerSec := %g\n", (double)NP * (double)steps_to_run / KernelSeconds);
    if(p.model == 0)
        printf("MaxAbsDeltaInvGammaMinusPz := %g\n", max_dgmpz);
    printf("Wrote %s (NP=%ld nT=%ld dTau=%g steps_run=%ld block=%d grid=%d)\n",
           out_path, NP, nT, dTau, steps_to_run, block_size, grid_size);

    return 0;
}
