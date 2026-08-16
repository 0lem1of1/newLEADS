#ifndef LEADS_STATE_H
#define LEADS_STATE_H

/*
  Binary particle-state dump shared between the Armadillo LEADS code and the
  standalone CPU/CUDA benchmark harness.

  Plain C++ only (no Armadillo, no GSL) so the same header compiles under both
  g++ and nvcc. All arrays are passed as raw double*; Armadillo mat(NP,3) is
  column-major, so posi.memptr() is already the [x0..xN-1, y0..yN-1, z0..zN-1]
  layout written here and no repacking is needed.

  Layout:
    LeadsStateHeader          72 bytes
    double TR   [nT]
    double posi [3*NP]        column-major
    double velo [3*NP]        column-major
    double momt [3*NP]        column-major
    double accl [3*NP]        column-major
    double egama[NP]

  File size = 72 + 8*(nT + 13*NP)
*/

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>

#define LEADS_STATE_MAGIC "LEADSST1"

struct LeadsStateHeader
{
    char    magic[8];      // "LEADSST1", not NUL-terminated
    int64_t NP;
    int64_t nT;
    double  dTau;
    double  t_shift;
    int64_t reserved[4];   // extend without a version bump
};

static_assert(sizeof(LeadsStateHeader) == 72, "LeadsStateHeader has unexpected padding");

struct LeadsState
{
    LeadsStateHeader hd;
    double *TR, *posi, *velo, *momt, *accl, *egama;
};

inline bool LeadsWriteState(const char* path,
                            int64_t NP, int64_t nT, double dTau, double t_shift,
                            const double* TR,
                            const double* posi, const double* velo,
                            const double* momt, const double* accl,
                            const double* egama)
{
    FILE* f = fopen(path,"wb");
    if(!f) return false;

    LeadsStateHeader hd;
    memset(&hd,0,sizeof(hd));
    memcpy(hd.magic,LEADS_STATE_MAGIC,8);
    hd.NP = NP; hd.nT = nT; hd.dTau = dTau; hd.t_shift = t_shift;

    const size_t n3 = (size_t)(3*NP);
    bool ok = fwrite(&hd,sizeof(hd),1,f) == 1;
    ok = ok && fwrite(TR   ,sizeof(double),(size_t)nT,f) == (size_t)nT;
    ok = ok && fwrite(posi ,sizeof(double),n3        ,f) == n3;
    ok = ok && fwrite(velo ,sizeof(double),n3        ,f) == n3;
    ok = ok && fwrite(momt ,sizeof(double),n3        ,f) == n3;
    ok = ok && fwrite(accl ,sizeof(double),n3        ,f) == n3;
    ok = ok && fwrite(egama,sizeof(double),(size_t)NP,f) == (size_t)NP;

    if(fclose(f) != 0) ok = false;
    return ok;
}

inline void LeadsFreeState(LeadsState& s)
{
    free(s.TR);   free(s.posi); free(s.velo);
    free(s.momt); free(s.accl); free(s.egama);
    s.TR = s.posi = s.velo = s.momt = s.accl = s.egama = NULL;
}

inline bool LeadsReadState(const char* path, LeadsState& s)
{
    memset(&s,0,sizeof(s));

    FILE* f = fopen(path,"rb");
    if(!f) return false;

    if(fread(&s.hd,sizeof(s.hd),1,f) != 1)            { fclose(f); return false; }
    if(memcmp(s.hd.magic,LEADS_STATE_MAGIC,8) != 0)    { fclose(f); return false; }
    if(s.hd.NP <= 0 || s.hd.nT <= 0)                    { fclose(f); return false; }

    const size_t nT = (size_t)s.hd.nT;
    const size_t NP = (size_t)s.hd.NP;
    const size_t n3 = 3*NP;

    s.TR    = (double*)malloc(nT*sizeof(double));
    s.posi  = (double*)malloc(n3*sizeof(double));
    s.velo  = (double*)malloc(n3*sizeof(double));
    s.momt  = (double*)malloc(n3*sizeof(double));
    s.accl  = (double*)malloc(n3*sizeof(double));
    s.egama = (double*)malloc(NP*sizeof(double));

    bool ok = s.TR && s.posi && s.velo && s.momt && s.accl && s.egama;
    ok = ok && fread(s.TR   ,sizeof(double),nT,f) == nT;
    ok = ok && fread(s.posi ,sizeof(double),n3,f) == n3;
    ok = ok && fread(s.velo ,sizeof(double),n3,f) == n3;
    ok = ok && fread(s.momt ,sizeof(double),n3,f) == n3;
    ok = ok && fread(s.accl ,sizeof(double),n3,f) == n3;
    ok = ok && fread(s.egama,sizeof(double),NP,f) == NP;

    fclose(f);
    if(!ok) LeadsFreeState(s);
    return ok;
}

#endif
