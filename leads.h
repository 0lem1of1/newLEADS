#include "leads_extra.h"
 
class LEADS_Input
{
public:
        LEADS_Input();
        //Simulation parameters
        bool RadiationReaction ; 
        int RRModel;
        double SimulationTime,TimeStep; 
        
        //LaserPulse
        double Intensity,Tau_FWHM,Wavelength,k,Omega,BeamWaist,a0;
        double PulseProfile,PulseModel,Polarization,PulsePhase;
        
        //ParticleBeam
        double ChargeParticle,MassParticle,GamaParticle;
        double PropTheta,PropPhi;
        int NoBeamPart;
        double GamaFWHM,BunchLength,BunchRadius,EnergyChirp;
        double x0,y0,z0;
 
        //Diagnostics
        double InterpolateLevel;
        double WriteTrajecAfter;
        bool CalculatePower,CalculateSpectrum,SignalEmitted,AngularDistribution;
        double DetecThetaStart,DetecThetaEnd,DetecDTheta;
        double DetecPhiStart,DetecPhiEnd,DetecDPhi;
        bool UseLogScale;
        double SpectrumStart,SpectrumEnd,SpectrumDO;
        
        bool WritePhaseSpace;int PhaseStoreInterval;
        bool WriteTrajectories;int EveryParticle;
        
        bool WriteSpatialFieldsProfile,RunFullSimulation; 
        double at_time;
        int xy,xz,yz;
         
        void ParameterInfo();
        double ChirpFunction(double,int); 

};

// Defining Laser Puslses
class LEADS_Laser
{
private:
	double Tau_FWHM,Omega;  
	double I0;	  
	int pola;	 
	double a0;       
	double w0;
	double delta,phase;
	int model;
	double profile;
	double zr,eps,k ;
	double DL_Const,E_Const,M_Const,L_Const;
public:
	LEADS_Input in;
	LEADS_Laser();
	void LaserProfile(double,rowvec,rowvec&,rowvec&);
	void Radial_Laser(double,rowvec,rowvec&,rowvec&);
	void Paraxial_Laser(double,rowvec,rowvec&,rowvec&);
	void CSPSW_Laser(double,rowvec,rowvec&,rowvec&);
	void PlainWave_Laser(double,rowvec,rowvec&,rowvec&);

	double w(double z){return w0 *sqrt(1 + SQ(z/zr));}
	double R(double z){return z + SQ(zr)/z ;}
	double Sn(int n,double z,double Psi,double PsiG)
	      {return PO(w0/w(z),n*1.0) * sin(Psi + n * PsiG);}
	double Cn(int n,double z,double Psi,double PsiG)
	      {return PO(w0/w(z),n*1.0) * cos(Psi + n * PsiG);}
};

// Dynamics
class LEADS_Dynamics 
{
 private: 
  	int NP,nT;
        double q_part,m_part;
        double parTheta,parPhi,t_shift,k,Omega;
        double RRFLAG,dTau,Time,SimuTime;
        double x_old,y_old,z_old;
        double RR_CONST_LL,RR_CONST_SOKOLOV,RR_CONST_FORD;
  	
  	mat posi,velo,momt,accl,bfield,efield,vT,vKAPPA;
  	cube vRP,vBETA,vBETADOT,vNNU;
  	vec TR,egama;
  	
  	double th_st,th_en,dTh,ph_st,ph_en,dPh;
  	vec dPdOmega; 
  	mat E_Signal,B_Signal;
  	mat AngDistri;
  	
 public:
   	int Num_Threads;
   	
   	LEADS_Input in;
 	LEADS_Laser la;
 	LEADS_Dynamics();
  	
 	void Run();
 	void InitializeParticles();
 	void WriteLaserPulse();
        void CalForceMovePart(double,int);
  	void CalRadiation();
        void RadiationModule(double,double);
        void SpectrumCalculate(double,double);
        void CalcEmittedSignal(double,double);
        double TR_Integrate(vec);
        Complex SP_Integrate(cx_vec);
        
        void WriteSpatialFields();
        void WriteSelectedData();
        void WritePhaseSpace();
        void WriteBunchTrajectories();
        void WritePlotInfo(string,int,int);   
};


class LEADS_Simulation
{
public:
      LEADS_Input in;
      LEADS_Dynamics dy;
 
      void Greetings();
      void Requirements();
      void ParameterChecks();
      void Run_Simulation(int);

};
