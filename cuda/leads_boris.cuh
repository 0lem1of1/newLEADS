#ifndef LEADS_BORIS_CUH
#define LEADS_BORIS_CUH

/*
  Portable core of the raw (no radiation-reaction) Boris pusher and the
  Model-0 ("Plain Wave") laser field, ported from:

    leads_traj.cpp  :: CalForceMovePart, final `else` branch (RRFLAG == 0)
    leads_laser.cpp :: LaserProfile, case 0

  Plain scalar C++ only (no Armadillo, no GSL) so this single header compiles
  unchanged under g++ (host-only build, cuda/bench_cpu.cpp) and nvcc
  (host+device build, cuda/bench_gpu.cu). LEADS_HD expands to
  `__host__ __device__` under nvcc and to nothing under g++.

  Every arithmetic expression below is written to match the exact operation
  order Armadillo's rowvec/cdot expressions produce, not just the same
  mathematical result — this is required for the CPU/GPU ports to reproduce
  the Armadillo reference to near machine precision (see cuda/statediff.cpp).
  In particular:
    - dot products are folded left-to-right as ((a*a)+(b*b))+(c*c), matching
      Armadillo's cdot() accumulation order for a 3-element rowvec.
    - `UP = term1*UM + term2*UxB + term2*del*udb*B` is grouped as
      `(term1*UM + term2*UxB) + (term2*del*udb*B)`, with the scalar prefactor
      on the last term built left-to-right as `((term2*del)*udb)*B` before it
      ever touches a vector component.
    - Model 0 does NOT divide by E_Const/M_Const — that normalization only
      applies to models 1 and 3 in leads_laser.cpp, so a templated
      laser_profile<MODEL>() must not carry it into the MODEL==0 case.
*/

#include <cmath>
#include "../leads_params.h"

#ifdef __CUDACC__
#define LEADS_HD __host__ __device__
#else
#define LEADS_HD
#endif

// Primary template intentionally left undefined: only MODEL==0 (Plain Wave)
// is ported. Instantiating any other MODEL is a link-time error by design.
template<int MODEL>
LEADS_HD void laser_profile(double t, double x, double y, double z,
                             const LeadsParams &p,
                             double &Ex, double &Ey, double &Ez,
                             double &Bx, double &By, double &Bz);

template<>
LEADS_HD inline void laser_profile<0>(double t, double /*x*/, double /*y*/, double z,
                                       const LeadsParams &p,
                                       double &Ex, double &Ey, double &Ez,
                                       double &Bx, double &By, double &Bz)
{
    // leads_laser.cpp case 0 has no x,y (transverse) dependence at all.
    const double eta = t - z + p.phase;

    // Chirping function is unused (hardcoded 0) in leads_laser.cpp's case 0.
    const double feta = 0.0, d_feta = 0.0;

    const double E0 = p.a0;
    const double B0 = p.a0;

    const double order = p.order;
    const double alpha = p.alpha;
    const double g_eta      = exp(-alpha * pow(fabs(eta), order));
    const double derivative = alpha * order * pow(fabs(eta), order - 1);

    const double delta   = p.delta;
    const double Deta_Dt = 1.0;
    const double Deta_Dz = -1.0;

    Ex = E0 * delta * g_eta * Deta_Dt
            * ( (1 + d_feta) * sin(eta + feta) + derivative * cos(eta + feta) );
    Ey = E0 * sqrt(1 - delta*delta) * g_eta * Deta_Dt
            * ( derivative*sin(eta+feta) - (1 + d_feta)*cos(eta+feta) );
    Ez = 0.0;

    Bx = B0 * sqrt(1 - delta*delta) * g_eta * Deta_Dz
            * ( derivative * sin(eta+feta) - (1 + d_feta)*cos(eta+feta) );
    By = -B0 * delta * g_eta * Deta_Dz
             * ( (1 + d_feta)*sin(eta+feta) + derivative * cos(eta+feta) );
    Bz = 0.0;
}

// Raw Boris push (radiation reaction off), transcribed scalar-for-scalar from
// leads_traj.cpp::CalForceMovePart's final `else` branch, including the
// magnetic-rotation preamble shared by every RRFLAG branch in that function.
//
// x,y,z / px,py,pz / vx,vy,vz / ax,ay,az are read as the particle's state at
// the START of this step and overwritten with its state at the END of it.
// vx,vy,vz on input is Vold (needed only to compute the output ax,ay,az
// finite-difference diagnostic, exactly as leads_traj.cpp does).
LEADS_HD inline void boris_step(
    double dTau, double q_part, double m_part,
    double &x, double &y, double &z,
    double &px, double &py, double &pz,
    double &vx, double &vy, double &vz,
    double &ax, double &ay, double &az,
    double Ex, double Ey, double Ez,
    double Bx, double By, double Bz)
{
    const double Vx_old = vx, Vy_old = vy, Vz_old = vz;

    double Px = px, Py = py, Pz = pz;

    double gamma = sqrt(1 + (((Px*Px) + (Py*Py)) + (Pz*Pz)));

    const double del_s = (q_part/m_part) * dTau/2.0;
    const double del   = del_s/gamma;
    const double del2  = del*del;
    const double B2    = ((Bx*Bx) + (By*By)) + (Bz*Bz);
    const double term1 = (1 - del2*B2) / (1 + del2*B2);
    const double term2 = 2*del / (1 + del2*B2);

    // UM = P + del_s * E
    const double UMx = Px + del_s*Ex;
    const double UMy = Py + del_s*Ey;
    const double UMz = Pz + del_s*Ez;

    // UxB = Cross(UM, B), matching leads_extra.h's Cross()
    const double UxBx = UMy*Bz - UMz*By;
    const double UxBy = UMz*Bx - UMx*Bz;
    const double UxBz = UMx*By - UMy*Bx;

    const double udb = ((UMx*Bx) + (UMy*By)) + (UMz*Bz);

    // UP = term1*UM + term2*UxB + term2*del*udb*B
    const double coeff = (term2*del)*udb;
    const double UPx = (term1*UMx + term2*UxBx) + coeff*Bx;
    const double UPy = (term1*UMy + term2*UxBy) + coeff*By;
    const double UPz = (term1*UMz + term2*UxBz) + coeff*Bz;

    // P = UP + del_s * E
    Px = UPx + del_s*Ex;
    Py = UPy + del_s*Ey;
    Pz = UPz + del_s*Ez;

    // Raw pusher (RRFLAG == 0): pbeta = P, no radiation-reaction correction.
    px = Px; py = Py; pz = Pz;

    gamma = sqrt(1 + (((px*px) + (py*py)) + (pz*pz)));

    vx = px/gamma; vy = py/gamma; vz = pz/gamma;

    x += vx*dTau; y += vy*dTau; z += vz*dTau;

    ax = (vx - Vx_old)/dTau;
    ay = (vy - Vy_old)/dTau;
    az = (vz - Vz_old)/dTau;
}

#endif
