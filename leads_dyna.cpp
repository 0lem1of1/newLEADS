#include "leads.h"

LEADS_Dynamics ::LEADS_Dynamics()
{
 NP       = in.NoBeamPart;
 parTheta = in.PropTheta * pi/180.0;
 parPhi   = in.PropPhi * pi/180.0;
 k        = 2 *pi /in.Wavelength;
 Omega    = in.Omega;
 dTau     = in.TimeStep;
 q_part   = in.ChargeParticle;
 m_part   = in.MassParticle;
 
 RRFLAG   = in.RadiationReaction > 0 ? 1 : 0; 
 RR_CONST_LL =   2 * k * re/3.0 ;
 RR_CONST_SOKOLOV =  Omega * tau0;
 RR_CONST_FORD = Omega * tau0;
 
 InitializeParticles();
  
 nT = int(2 * t_shift / dTau) + 1;
 if(in.StoreTrajectories)
 {
  vRP   = zeros(NP,3,nT);
  vBETA = zeros(NP,3,nT);
  vBETADOT = zeros(NP,3,nT);
  vNNU = zeros(NP,3,nT);
  vT = zeros(NP,nT);
  vKAPPA = zeros(NP,nT);
 }
 traj0 = zeros(nT,3);
 beta0 = zeros(nT,3);
 TR = linspace(-t_shift,t_shift,nT);
  
 th_st = in.DetecThetaStart * pi/180.0;
 th_en = in.DetecThetaEnd   * pi/180.0; 
 dTh   = in.DetecDTheta     * pi/180.0;
 ph_st = in.DetecPhiStart   * pi/180.0;
 ph_en = in.DetecPhiEnd     * pi/180.0;
 dPh   = in.DetecDPhi       * pi/180.0;
  
}

void LEADS_Dynamics ::Run()
{
  WriteLaserPulse();
  LEADS_Dynamics();
  
  Print("Starting interaction dynamics");
  auto PusherStart = high_resolution_clock::now();
  for(int i = 0; i < nT ; i++)
  {
    for(int j = 0;j< NP;j++)
    {
     CalForceMovePart(TR(i),j);
    }
    traj0.row(i) = posi.row(0);
    beta0.row(i) = velo.row(0);
    if(in.StoreTrajectories)
    {
     vRP.slice(i) = posi;
     vBETA.slice(i) = velo;
     vBETADOT.slice(i) = accl;
    }
  }
  auto PusherEnd = high_resolution_clock::now();
  double PusherSeconds = duration_cast<duration<double>>(PusherEnd - PusherStart).count();
  VPrint("PusherLoopTime_sec",PusherSeconds);
  VPrint("FinalStateChecksum",accu(abs(posi)) + accu(abs(momt)));
  Print("Elec. trajectories calculated");
  
  WriteSelectedData();
  
  CalRadiation();
}


