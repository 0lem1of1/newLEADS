#include "leads.h"

void LEADS_Simulation ::  ParameterChecks()
{
 Print(in.RadiationReaction > 0 ? "Radiation Reaction Included" : "Radiation Reaction Excluded") ;
 
 if(in.RadiationReaction > 0 && in.RRModel <= 3)
 {
  switch(in.RRModel)
  {
   case 0:VPrint("RR Model","LL without derivative");break;
   case 1:VPrint("RR Model","LL with derivative");break; 
   case 2:VPrint("RR Model","Sokolov");break; 
   case 3:VPrint("RR Model","Ford & Connell");break; 
  }
 }
 
 if(in.RadiationReaction > 0 && in.RRModel > 3)
 {
  Print("Invalid RRModel, aborting LEADS");cout << endl;
  exit(0);
 }
  
 
 if(in.PulseModel > 3 )
 {
  Print("Invalid PulseModel, aborting LEADS");cout << endl;
  exit(0);
 }
 else
 {
   switch(int(in.PulseModel))
   {
    case 0:VPrint("Laser Model","Plain Wave");break;
    case 1:VPrint("Laser Model","Paraxial Beam");break; 
    case 2:VPrint("Laser Model","CSPSW Beam");break; 
    case 3:VPrint("Laser Model","Radially Polarized");break; 
   }
 }
 
 if(in.Polarization > 1 )
 {
  Print("Invalid Polarization, aborting LEADS");cout << endl;
  exit(0);
 }

 if(!in.StoreTrajectories && (in.CalculateSpectrum || in.AngularDistribution ||
    in.SignalEmitted || in.CalculatePower || in.WritePhaseSpace || in.WriteTrajectories))
 {
  Print("StoreTrajectories=false requires all trajectory diagnostics off, aborting LEADS");cout << endl;
  exit(0);
 }

}

void LEADS_Simulation :: Run_Simulation(int NThreads)
{
    Requirements();
    
    Print("Doing parameter checks");
    ParameterChecks();
    
    Print("Starting simulation");
    
    dy.Num_Threads = NThreads;
    dy.Run();	 
}
