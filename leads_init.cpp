#include "leads.h"

void LEADS_Dynamics ::InitializeParticles()
{
 Print("Initializing electron bunch");
 
 posi = zeros(NP,3);    // position x,y,z
 velo = zeros(NP,3);    // velocity vx,vy,vz
 momt = zeros(NP,3);    // momentum px,py,pz
 accl = zeros(NP,3);    // acceleation ax,ay,az
 bfield = zeros(NP,3);  // magentic field:needed for RR 
 efield = zeros(NP,3);  // electric field:needed for RR 
 egama = zeros(NP);     // gama of particles
 
 const gsl_rng_type * temp;
 temp = gsl_rng_default;

 gsl_rng *Gama, *Length, *Radius; 
 gsl_rng_env_setup();
 Gama   = gsl_rng_alloc (temp);
 Length = gsl_rng_alloc (temp);
 Radius = gsl_rng_alloc (temp);
 
 double E_Chirp  = -in.EnergyChirp;
 double ch_temp0 = 2.0 * abs(E_Chirp) / in.BunchLength; 
 double muGama = in.GamaFWHM;
 double muLeng = in.BunchLength;
 double muRadi = in.BunchRadius;
 double x_temp,y_temp,z_temp,vel_temp;
 int FRx,FRy,FRz; // flags for radius along x,y,z
 int FLx,FLy,FLz; // flags for length along x,y,z
 
 if(int(in.PropTheta) == 180 && int(in.PropPhi) == 0)
 {FLz = 1; FLx = FLy = FRz = 0; FRx = FRy = 1;  }

 if(abs(int(in.PropTheta)) == 90 && int(in.PropPhi) == 0)
 {FLx = 1; FLz = FLy = FRx = 0; FRz = FRy = 1; }

 if(abs(int(in.PropTheta)) == 90 && abs(int(in.PropPhi)) == 90)
 {FLy = 1; FLx = FLz = FRy = 0; FRx = FRz = 1; }
 
  
 for(int i = 0;i< NP;i++)
 {
  if(i == 0)
  {
   posi(i,0) = in.z0 * sin(pi+parTheta)*cos(parPhi) + in.x0;
   posi(i,1) = in.z0 * sin(pi+parTheta)*sin(parPhi) + in.y0;
   posi(i,2) = in.z0 * cos(pi+parTheta);
   egama(i) = in.GamaParticle;
  }
  else
  {
    x_temp = in.x0 + FRx * gsl_ran_gaussian (Radius, muRadi) 
                   + FLx * (gsl_ran_flat (Length,in.x0+0.5*muLeng, in.x0-0.5*muLeng) - in.x0);
    y_temp = in.y0 + FRy * gsl_ran_gaussian (Radius, muRadi) 
                   + FLy * (gsl_ran_flat (Length,-0.5*muLeng,0.5*muLeng) );
    z_temp = in.z0 + FRz * gsl_ran_gaussian (Radius, muRadi) 
                   + FLz * (gsl_ran_flat (Length,in.z0+0.5*muLeng, in.z0-0.5*muLeng) - in.z0);

    posi(i,0) = x_temp + z_temp * sin(pi+parTheta)*cos(parPhi);
    posi(i,1) = y_temp + z_temp * sin(pi+parTheta)*sin(parPhi);
    posi(i,2) = z_temp * cos(pi+parTheta) ;  
    
    // energy chirp of electron bunch, only along z direction propagation
   egama(i) = in.GamaParticle * (1.0 + FLz * Sign(E_Chirp) * ch_temp0 * (-in.z0 + posi(i,2)) );
  }
  
  vel_temp = sqrt(1.0 - 1/SQ(egama(i)));
  velo(i,0) = vel_temp * sin(parTheta)*cos(parPhi);
  velo(i,1) = vel_temp * sin(parTheta)*sin(parPhi);
  velo(i,2) = vel_temp * cos(parTheta);
  momt(i,0) = velo(i,0) * egama(i); 
  momt(i,1) = velo(i,1) * egama(i); 
  momt(i,2) = velo(i,2) * egama(i);  
 }
 
 //posi.clean(1.0E-13);
 //velo.clean(1.0E-13);
 //momt.clean(1.0E-13);
 
 double t_shift0 = sqrt(cdot(posi.row(0),posi.row(0)))/sqrt(cdot(velo.row(0),velo.row(0)));
 double t_shift1 = 0.5 * in.SimulationTime * in.Tau_FWHM;
 
 t_shift = in.GamaParticle > 1 ? t_shift0 : t_shift1;

  
 Print("Saving initial electron bunch"); 
 ofstream ebunch("./Data/eBunch_t0.dat",ios::out);
 ebunch << "#x,y,z,gama" << endl;
 mat store = posi;
 store.insert_cols(3,egama);
 store.save(ebunch,raw_ascii);
 store.clear(); 
 ebunch.close();
  
 WritePlotInfo("ebunch",0,0); 
}
