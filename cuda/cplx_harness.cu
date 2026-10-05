/*
  Host-vs-device agreement harness for the Model 2 (focused 3D pulse) field.

  Evaluates the exact expression sequence from leads_laser.cpp case 2 three
  ways at the same sampled points:

    A  host,   std::complex<double>      (what Armadillo-side LEADS uses)
    B  host,   thrust::complex<double>
    C  device, thrust::complex<double>   (what bench_gpu would use)

  and reports field-level (Ex,Ey,Ez,Bx,By,Bz) disagreement A-vs-B (library
  difference only) and A-vs-C (library + device codegen), plus how often the
  `Rc = imag(Rc)<0 ? -Rc : Rc` branch selection differs.

  Error metric: |a-b| / s, where s = max over the 6 components of |A| at that
  point (point-local scale, so fields crossing zero don't blow up the ratio).
  Points where all six reference components are exactly 0 are skipped.
  The point-local metric is dominated by cancellation noise in the far tails
  (fields ~1e-13 of the peak), where even std::complex is off vs a
  long-double reference by the same amount -- so the pass/fail number is
  max|abs error| / global max|field| (printed as max|abs|/fieldmax).

  Usage: cplx_harness [--params file] [--pola N] [--a0 X] [--a0 X] [--n N] [--seed S]
  Build: see cuda/Makefile (cplx_harness / cplx_harness_nofma).
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <random>
#include <thrust/complex.h>
#include "../leads_params.h"

#define HD __host__ __device__

// ---- tiny adapter layer so one template serves both complex types ---------
template<class T> HD inline std::complex<T>    csqrt(const std::complex<T>& z)    { return std::sqrt(z); }
template<class T> HD inline std::complex<T>    cexp (const std::complex<T>& z)    { return std::exp(z); }
template<class T> HD inline thrust::complex<T> csqrt(const thrust::complex<T>& z) { return thrust::sqrt(z); }
template<class T> HD inline thrust::complex<T> cexp (const thrust::complex<T>& z) { return thrust::exp(z); }

struct Out { double Ex,Ey,Ez,Bx,By,Bz; int rc_flipped; double rcre, rcim; };

// Verbatim transcription of leads_laser.cpp case 2 (operation order kept).
// R is the real scalar type of Complex; every literal/scalar goes through R
// so the same body also instantiates with std::complex<long double>, which
// serves as a higher-precision reference to judge who is right where the
// double-precision paths disagree.
template<class Complex, class R = typename Complex::value_type>
HD inline Out model2(double t_, double x_, double y_, double z_, const LeadsParams& p)
{
    const R t = t_, x = x_, y = y_, z = z_;
    const R one = 1, two = 2, half = 0.5;
    const Complex iota0(0,1);
    const R E0 = p.a0;
    R z0 = R(p.k)*R(p.zr);
    R t0 = p.phase;
    R zeta2 = (R)p.pola;
    R T = R(p.Tau_FWHM)/sqrt(R(8)*log(R(2)));
    Complex Rc = csqrt(Complex(x*x + y*y,0) + (z + iota0*z0)*(z + iota0*z0));
    int flipped = 0;
    if(Rc.imag() < 0) { Rc = -Rc; flipped = 1; }
    Complex tc = (t - t0 + iota0*z0);
    Complex tau = (tc - Rc);
    Complex p0 = (z0 * E0)/sqrt((one - one/z0 + one/(T*T) + one/(z0*z0)) * (one - one/z0 + one/(T*T)));
    R phi0 = 0;
    Complex phase0 = cexp(iota0*(tau + phi0));
    Complex P_0 = p0*cexp(-half*(tau*tau/T/T)) * phase0;
    Complex Zc = z + iota0*z0;

    Complex f = (one + (iota0*tau)/(T*T))*(one + (iota0*tau)/(T*T))
                - (one/(Rc*Rc))*(one - tc*Rc/(T*T) + iota0*Rc);
    Complex g = -f + (two/(Rc*Rc))*(one - tau*Rc/(T*T) + iota0*Rc);
    Complex hf = f + one/(Rc*Rc);

    Complex dEx = (P_0)/(Rc*sqrt(one+zeta2*zeta2))*(f + g*x*(x + iota0*zeta2*y)/(Rc*Rc));
    Complex dEy = (P_0)/(Rc*sqrt(one+zeta2*zeta2))*(iota0*f*zeta2 + g*y*(x + iota0*zeta2*y)/(Rc*Rc));
    Complex dEz = (P_0)/(Rc*sqrt(one+zeta2*zeta2))*(g*Zc*(x + iota0*zeta2*y)/(Rc*Rc));

    R factor = 1;
    Complex dBx = (P_0*hf)/(Rc*Rc*sqrt(one+zeta2*zeta2)*factor)*(-iota0*zeta2*Zc);
    Complex dBy = (P_0*hf)/(Rc*Rc*sqrt(one+zeta2*zeta2)*factor)*(Zc);
    Complex dBz = (P_0*hf)/(Rc*Rc*sqrt(one+zeta2*zeta2)*factor)*(iota0*zeta2*x - y);

    Out o;
    o.Ex = (double)dEx.real(); o.Ey = (double)dEy.real(); o.Ez = (double)dEz.real();
    o.Bx = (double)dBx.real(); o.By = (double)dBy.real(); o.Bz = (double)dBz.real();
    o.rc_flipped = flipped; o.rcre = (double)Rc.real(); o.rcim = (double)Rc.imag();
    return o;
}

struct Pt { double t,x,y,z; };

__global__ void kern(const Pt* pts, Out* out, int n, LeadsParams p)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i < n) out[i] = model2<thrust::complex<double>>(pts[i].t,pts[i].x,pts[i].y,pts[i].z,p);
}

#define CK(x) do{ cudaError_t e_=(x); if(e_!=cudaSuccess){ fprintf(stderr,"CUDA: %s (%s:%d)\n",cudaGetErrorString(e_),__FILE__,__LINE__); exit(2);} }while(0)

static inline void comps(const Out& o, double c[6]) { c[0]=o.Ex;c[1]=o.Ey;c[2]=o.Ez;c[3]=o.Bx;c[4]=o.By;c[5]=o.Bz; }

static double G_FMAX = 1.0;
struct Stats { double maxabs=0, maxe=0, sum=0; long n=0, gross=0, nonfinite=0, rcdiff=0; std::vector<double> e; };

static void compare(const char* name, const std::vector<Out>& A, const std::vector<Out>& B,
                    const std::vector<Pt>& pts, const char* tag)
{
    Stats s;
    s.e.reserve(A.size());
    long worst = -1;
    for(size_t i=0;i<A.size();i++)
    {
        double a[6], b[6]; comps(A[i],a); comps(B[i],b);
        double sc = 0; bool fin = true;
        for(int j=0;j<6;j++){ sc = std::max(sc, fabs(a[j])); if(!std::isfinite(a[j])||!std::isfinite(b[j])) fin=false; }
        if(!fin){ s.nonfinite++; continue; }
        if(A[i].rc_flipped != B[i].rc_flipped) s.rcdiff++;
        if(sc == 0) continue;
        double m = 0;
        for(int j=0;j<6;j++){ m = std::max(m, fabs(a[j]-b[j])/sc); s.maxabs = std::max(s.maxabs, fabs(a[j]-b[j])); }
        s.e.push_back(m); s.sum += m; s.n++;
        if(m > s.maxe){ s.maxe = m; worst = (long)i; }
        if(m > 1e-6) s.gross++;
    }
    std::sort(s.e.begin(), s.e.end());
    auto q = [&](double f){ return s.e.empty()?0.0:s.e[std::min(s.e.size()-1,(size_t)(f*s.e.size()))]; };
    printf("%-28s [%s] n=%ld  mean=%.3e  p50=%.3e  p99=%.3e  p99.9=%.3e  max=%.3e  max|abs|/fieldmax=%.3e  >1e-6: %ld  nonfinite: %ld  Rc-branch-diff: %ld\n",
           name, tag, s.n, s.n?s.sum/s.n:0.0, q(0.5), q(0.99), q(0.999), s.maxe, s.maxabs/G_FMAX, s.gross, s.nonfinite, s.rcdiff);
    if(worst >= 0)
    {
        const Pt& w = pts[worst];
        printf("    worst point: t=%.6g x=%.6g y=%.6g z=%.6g\n", w.t,w.x,w.y,w.z);
    }
}

int main(int argc, char** argv)
{
    const char* ppath = "../TestCase/bench/params.bin";
    int  pola_override = -999;
    double a0_override = -1;
    long N = 1000000;
    unsigned seed = 12345;
    for(int i=1;i<argc;i++)
    {
        if(!strcmp(argv[i],"--params") && i+1<argc) ppath = argv[++i];
        else if(!strcmp(argv[i],"--pola") && i+1<argc) pola_override = atoi(argv[++i]);
        else if(!strcmp(argv[i],"--a0") && i+1<argc) a0_override = atof(argv[++i]);
        else if(!strcmp(argv[i],"--n") && i+1<argc) N = atol(argv[++i]);
        else if(!strcmp(argv[i],"--seed") && i+1<argc) seed = (unsigned)atol(argv[++i]);
        else { fprintf(stderr,"usage: %s [--params f] [--pola N] [--a0 X] [--n N] [--seed S]\n",argv[0]); return 1; }
    }

    LeadsParams p;
    if(!LeadsReadParams(ppath,p)){ fprintf(stderr,"cannot read params %s\n",ppath); return 1; }
    if(pola_override != -999) p.pola = pola_override;
    if(a0_override >= 0) p.a0 = a0_override;

    const double z0 = p.k*p.zr, w0n = p.k*p.w0;
    const double T  = p.Tau_FWHM/sqrt(8*log(2.0));
    printf("params: k=%.6g w0(norm)=%.6g zr(norm)=z0=%.6g Tau_FWHM=%.6g T=%.6g a0=%.6g phase=%.6g pola=%d (file model=%d)\n",
           p.k, w0n, z0, p.Tau_FWHM, T, p.a0, p.phase, p.pola, p.model);

    // Sample points: pulse-relevant times, focal-region positions, plus
    // forced z=0 plane (20%) and on-axis x=y=0 (10%) subsets.
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(-1.0,1.0);
    std::normal_distribution<double> Nd(0.0,1.0);
    std::vector<Pt> pts(N);
    for(long i=0;i<N;i++)
    {
        Pt& q = pts[i];
        q.x = 3*w0n*U(rng); q.y = 3*w0n*U(rng); q.z = 3*z0*U(rng);
        int r = (int)(i % 10);
        if(r < 2) q.z = 0.0;           // z=0 plane
        if(r == 2) { q.x = 0; q.y = 0; } // on axis
        if(r == 3) { q.x = 0; q.y = 0; q.z = 0; } // origin
        double rr = sqrt(q.x*q.x + q.y*q.y + q.z*q.z);
        q.t = rr + p.phase + 3*T*Nd(rng); // near/inside the pulse envelope
    }

    std::vector<Out> A(N), B(N), C(N), L(N);
    for(long i=0;i<N;i++) A[i] = model2<std::complex<double>>(pts[i].t,pts[i].x,pts[i].y,pts[i].z,p);
    for(long i=0;i<N;i++) B[i] = model2<thrust::complex<double>>(pts[i].t,pts[i].x,pts[i].y,pts[i].z,p);

    for(long i=0;i<N;i++) L[i] = model2<std::complex<long double>>(pts[i].t,pts[i].x,pts[i].y,pts[i].z,p);

    Pt* dP; Out* dO;
    CK(cudaMalloc(&dP, N*sizeof(Pt))); CK(cudaMalloc(&dO, N*sizeof(Out)));
    CK(cudaMemcpy(dP, pts.data(), N*sizeof(Pt), cudaMemcpyHostToDevice));
    kern<<<(N+255)/256,256>>>(dP,dO,(int)N,p);
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(C.data(), dO, N*sizeof(Out), cudaMemcpyDeviceToHost));
    cudaFree(dP); cudaFree(dO);

    // field magnitude overview so error numbers can be read in context
    double fmax = 0; long zeros = 0;
    for(long i=0;i<N;i++){ double c[6]; comps(A[i],c); double s=0; for(int j=0;j<6;j++) s=std::max(s,fabs(c[j])); fmax=std::max(fmax,s); if(s==0) zeros++; }
    G_FMAX = fmax;
    printf("N=%ld  max|field|=%.4g  exactly-zero points=%ld\n\n", N, fmax, zeros);

    // who is right? each double path vs the long-double reference
    compare("std host vs LONGDOUBLE ref", L, A, pts, "all");
    compare("thrust host vs LONGDOUBLE", L, B, pts, "all");
    compare("thrust DEVICE vs LONGDOUBLE", L, C, pts, "all");
    printf("\n");

    // all points, then the risky subsets separately
    compare("std host vs thrust host", A, B, pts, "all");
    compare("std host vs thrust DEVICE", A, C, pts, "all");
    compare("thrust host vs DEVICE", B, C, pts, "all");

    for(int sub=0; sub<2; sub++)
    {
        std::vector<Out> a,b,c; std::vector<Pt> pp;
        for(long i=0;i<N;i++)
        {
            int r = (int)(i%10);
            bool in = sub==0 ? (r<2 || r==3) : (r==2 || r==3);
            if(in){ a.push_back(A[i]); b.push_back(B[i]); c.push_back(C[i]); pp.push_back(pts[i]); }
        }
        const char* tag = sub==0 ? "z=0 plane" : "on-axis";
        compare("std host vs thrust host", a, b, pp, tag);
        compare("std host vs thrust DEVICE", a, c, pp, tag);
    }
    return 0;
}
