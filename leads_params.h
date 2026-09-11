#ifndef LEADS_PARAMS_H
#define LEADS_PARAMS_H

/*
  Scalar-parameter sidecar shared between the Armadillo LEADS code and the
  standalone CPU/CUDA Boris-pusher benchmark harness (see leads_state.h for
  the companion particle-state binary dump).

  Plain C++ only (no Armadillo, no GSL) so the same header compiles under
  both g++ and nvcc. LEADS_Laser::FillParams() populates this from the same
  private members / inline expressions used by LaserProfile(), so the
  sidecar always matches what the field calculation actually used.

  Layout:
    char    magic[8]                "LEADSPR1", not NUL-terminated
    double  k,Omega,Tau_FWHM,a0,w0,delta,phase,zr,eps,profile
    double  DL_Const,E_Const,M_Const,L_Const
    double  order,alpha
    double  q_part,m_part,dTau
    int32_t model,pola
    int64_t reserved[4]              extend without a version bump
*/

#include <cstdio>
#include <cstdint>
#include <cstring>

#define LEADS_PARAMS_MAGIC "LEADSPR1"

struct LeadsParams
{
    char    magic[8];

    double  k, Omega, Tau_FWHM, a0, w0, delta, phase, zr, eps, profile;
    double  DL_Const, E_Const, M_Const, L_Const;
    double  order, alpha;
    double  q_part, m_part, dTau;

    int32_t model, pola;

    int64_t reserved[4];
};

static_assert(sizeof(LeadsParams) == 200, "LeadsParams has unexpected padding");

inline bool LeadsWriteParams(const char* path, const LeadsParams& p_in)
{
    LeadsParams p = p_in;
    memset(p.magic,0,sizeof(p.magic));
    memcpy(p.magic,LEADS_PARAMS_MAGIC,8);
    memset(p.reserved,0,sizeof(p.reserved));

    FILE* f = fopen(path,"wb");
    if(!f) return false;

    bool ok = fwrite(&p,sizeof(p),1,f) == 1;

    if(fclose(f) != 0) ok = false;
    return ok;
}

inline bool LeadsReadParams(const char* path, LeadsParams& p)
{
    memset(&p,0,sizeof(p));

    FILE* f = fopen(path,"rb");
    if(!f) return false;

    bool ok = fread(&p,sizeof(p),1,f) == 1;
    fclose(f);

    if(!ok) return false;
    if(memcmp(p.magic,LEADS_PARAMS_MAGIC,8) != 0) return false;

    return true;
}

#endif
