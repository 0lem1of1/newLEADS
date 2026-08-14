#include "./files/leads.h"

LEADS_Input :: LEADS_Input()
{
  //------------------------------------------------------------------------------------------ 
  //RunParameters
   RadiationReaction = false;  // Include Radiation Reaction true/false
   RRModel           = 0;      // 0 <= RRModel <= 3; See Parameter info at bottom
   SimulationTime    = 3;      // Simulation time (x Tau_FWHM)
   TimeStep          = 5.0E-3; // Simulation time step in 'femtoseconds' or 'dimensionless'
   
  //------------------------------------------------------------------------------------------ 
  // Laser Parameters
   Intensity  = 2.0E+18;                                           // Intensity in 'W/cm2'
   Wavelength = 800.0E-9;                                          // Wavelength in 'nanometer'
   Omega      = 2.0*pi*c/Wavelength;                               // Frequency 
   a0         = 10.0;                                              // Laser amplitude
   Tau_FWHM   = 4 * (2*pi);                                        // Duration in cycles
   //a0       = sqrt(Intensity*SQ(Wavelength/1e-6)/1.36627476e18); // Laser amplitude
   //Tau_FWHM = 10.0 * (Omega * fs);                               // Duration in fs
   
   PulseProfile = log(2.0)/log(2.0);    // log(n)/log(2) ; n =2 gaussian, n = 4 super-gaussian etc. 
   PulseModel   = 1;    // 0 <= PulseModel <= 3; See Parameter info at bottom
   BeamWaist    = 3;    // Beam Waist size in 'micrometers'
   Polarization = 0;    // 0 - Linear (Ex,By,kz) ; 1 - Circular (ExBx,EyBy,kz)
   PulsePhase   = 0;    // in radians
   
  //------------------------------------------------------------------------------------------ 
  // Particle Beam Parameters
   NoBeamPart     = 1;     // Number of particles in beam
   ChargeParticle = -1.0;  // Charge of particle in units of 'e'
   MassParticle   = 1.0;   // Mass of particle in units of 'me' 
   GamaParticle   = 500.0; // Initial gama of the particle

   PropTheta  = 180.0;    // Direction of propagation (Angle with z-axis)
   PropPhi    = 0.0;      // Direction of propagation (Angle with x-axis) 
	     
   GamaFWHM    = 1;       // FWHM in Gama for Gaussian distribution of particle beam
   BunchLength = 4;       // Length of cyclindrical bunch (will count from z0/x0)
   BunchRadius = 1;       // Radius of cylindrical bunch (will be centred around x0,y0) 
   //BunchRadius = BunchLength/3.0 ; // Radius of cylindrical bunch (will be centred around x0,y0) 
   
   EnergyChirp = 0 ;    // energy chirp in fraction of GamaParticle  
                              
   x0  = 0.0;   // inital position of the particle(s)
   y0  = 0.0;   // inital position of the particle(s)
   z0  = 20.0; // inital position of the particle(s)
  
  //------------------------------------------------------------------------------------------ 
  //Diagnostics: Spectrum Module
   WriteTrajecAfter  = 1;      // Write Phase space after time steps
   CalculatePower    = false;   // true/false
   CalculateSpectrum = false;   // true/false
   SignalEmitted     = false;   // true/false
   AngularDistribution = false; // true/false
   
   InterpolateLevel  = 1.0;     // insert points between two adjacent points while integration  
    
   DetecThetaStart   = 180.0;   // Direction of radiation detection (Angle with z-axis)
   DetecThetaEnd     = 180.0;   // Direction of radiation detection (Angle with z-axis)
   DetecDTheta       = 10.0;     // dTheta

   DetecPhiStart     = 0.0;     // Direction of radiation detection (Angle with x-axis)
   DetecPhiEnd       = 0.0;     // Direction of radiation detection (Angle with x-axis)
   DetecDPhi         = 10.0;     // dPhi
   
   UseLogScale       = true;    // true/false (Spectrum will be in power of 10)

   SpectrumStart     = 3;       // Spectrum start in terms of w/w0
   SpectrumEnd       = 8;       // Spectrum end in terms of w/w0
   SpectrumDO        = 0.1;    // Spectrum increment in terms of w/w0 
   
   //------------------------------------------------------------------------------------------ 
   //Diagnostics: Phase space and trajectory module
   
   WritePhaseSpace = false;    // true/false: store phase space
   PhaseStoreInterval = 10;   // store phase space after time intervals 
   
   WriteTrajectories = false;  // true/false: store bunch trajectory in time
   EveryParticle = 10;        // for every # particle
   
   //------------------------------------------------------------------------------------------ 
   //Diagnostics: Spatial laser profiles
   
   WriteSpatialFieldsProfile = true;
   RunFullSimulation = true;
   StoreTrajectories = false;   // false => skip NP*3*nT history cubes, final-state only
   at_time = 0;
   xy = 1; xz = 0; yz = 0; 
 }
   
 void LEADS_Input :: ParameterInfo()
 {
  #if(0)
  
   RRModel:
       0 - Landau-Lifshitz Model, without derivative term
   	   Radiation reaction effects on radiation pressure acceleration
     	   M Tamburini et. al., New J. Phys. 12, 123005,(2010)
     	   https://doi.org/10.1088/1367-2630/12/12/123005
        
       1 - Landau-Lifshitz Model, with EM field time derivative terms
   	   Radiation reaction effects on radiation pressure acceleration
     	   M Tamburini et. al., New J. Phys. 12, 123005,(2010)
     	   https://doi.org/10.1088/1367-2630/12/12/123005
     	   
       2 - Dynamics of emitting electrons in strong laser fields
     	   I V Sokolov et. al., Phys. Plasmas 16, 093115 (2009)
           https://doi.org/10.1063/1.3236748
           
       3 - Relativistic form of Radiation Reaction
           G W Ford and R F OConnell, Phys. Letters A 174, 182-184, (1993)
           https://doi.org/10.1016/0375-9601(93)90755-O
      
   
   PulseProfile: 
   		log(n)/log(2) 
                n = 1 - Infinite Plane Wave
   		    2 - exp(-t^2)
   		    3 - exp(-t^3)
   		    4 - exp(-t^4)
   		    2.5 - exp(-t^2.5)
               fractional power of 2 is also possible 
   
   PulseModel:
   		0 - Plane Wave
   		
   		1 - Electron scattering and acceleration by a tightly focused laser beam
  		    Yousef I. Salamin, Guido R. Mocken, and Christoph H. Keitel
  		    Phys. Rev. ST Accel. Beams 5, 101301 - Published 18 October 2002
                   https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.5.101301
                    
   		2 - Subcycle Pulsed Focused Vector Beams
  		    Qiang Lin, Jian Zheng, and Wilhelm Becker
  		    Phys. Rev. Lett. 97, 253902 - Published 19 December 2006
                   https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.253902
  
               3 - Acceleration of proton bunches by petawatt chirped radially polarized laser pulses
 		   Jian-Xing Li, Yousef I. Salamin, Benjamin J. Galow, and Christoph H. Keitel
                   Phys. Rev. A 85, 063832 - Published 25 June 2012
 		   https://journals.aps.org/pra/abstract/10.1103/PhysRevA.85.063832
 #endif
 
 }
