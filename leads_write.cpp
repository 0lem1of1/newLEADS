#include "leads.h"

void LEADS_Dynamics ::WritePlotInfo(string pinfo,int NCols,int NRows)
{
 
 string filename = "./Data/PlottingDetails.txt";
    
 if(pinfo == "signal")
 {
    ofstream opl(filename.c_str(),ios::app);
    opl << "\n Emitted signal for ang. range : " ; 
    opl << "\n File: '/Data/emitted_signal.dat' : ASCII";
    opl << "\n Columns: theta,phi,t_ret,t_lab,Ex,Ey,Ez,dP/dOmega";
    opl << "\n gnuplot: standard ASCII data plotting";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }
 
 if(pinfo == "laser")
 {
    ofstream opl(filename.c_str(),ios::app);
    opl << "\n Laser Pulse : " ; 
    opl << "\n File: '/Data/laser_pulse.dat' : ASCII";
    opl << "\n Columns: t,Ex,Ey,Ez";
    opl << "\n gnuplot: standard ASCII data plotting";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }
 
 if(pinfo == "ebunch")
 {
    ofstream opl(filename.c_str(),ios::app);
    opl << "\n Electron bunch at t = 0: " ; 
    opl << "\n File: '/Data/eBunch_t0.dat' : ASCII";
    opl << "\n Columns: x,y,z,gama";
    opl << "\n gnuplot: standard ASCII data plotting";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }  
 
 if(pinfo == "etraj")
 {
    ofstream opl(filename.c_str(),ios::app);
    opl << "\n Electron trajectory : " ; 
    opl << "\n File: '/Data/eTraj_n0.dat' : ASCII";
    opl << "\n Columns: t,x,y,z,vx,vy,vz,gama";
    opl << "\n gnuplot: standard ASCII data plotting";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 } 
 
 if(pinfo == "phase")
 {  
   ofstream opl(filename.c_str(),ios::app);
    opl << "\n Phase-Space plots : " ; 
    opl << "\n File: '/Data/phasespace/time-nxx' : BINARY";
    opl << "\n Columns: x,y,z,vx,vy,vz";
    opl << "\n gnuplot: following syntax need to be used";
    opl << "\n gnuplot> p 'time-nxx' binary format='%" << NCols << "lf' u 1:2 w lines";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }
 
 if(pinfo == "sfields")
 {  
   ofstream opl(filename.c_str(),ios::app);
    opl << "\n Spatial profiles of the fields : " ; 
    opl << "\n File: '/Data/fields/xx.dat' : MATRIX";
    opl << "\n gnuplot: following syntax need to be used";
    opl << "\n gnuplot> p 'xx.dat' matrix w image";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }
 
 if(pinfo == "butraj")
 {  
    ofstream opl(filename.c_str(),ios::app);
    opl << "\n Particle trajectories : " ; 
    opl << "\n File: '/Data/particles/part-xx' : BINARY";
    opl << "\n Columns: t,x,y,z,vx,vy,vz";
    opl << "\n gnuplot: following syntax need to be used";
    opl << "\n gnuplot> p 'part-xx' binary format='%" << NCols << "lf' u 1:2 w lines";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }
 
 if(pinfo == "spec")
 {  
    ofstream opl(filename.c_str(),ios::app);
    opl << "\n Radiation spectrum ang. range: " ; 
    opl << "\n File: '/Data/spectrum.dat' : ASCII";
    opl << "\n Columns: Theta, Phi, w/w0, d2I/dwdO";
    opl << "\n gnuplot: standard ASCII data plotting";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }  
 
 if(pinfo == "angdist")
 {  
    ofstream opl(filename.c_str(),ios::app);
    opl << "\n Angular distribution of energy radiated : " ; 
    opl << "\n File: '/Data/angl_distri.bin' : BINARY";
    opl << "\n Columns: matrix(theta,phi)";
    opl << "\n gnuplot: need following syntax";
    opl << "\n gnuplot> p 'angl_distri.bin' binary array=(" << NCols << "," << NRows 
                 << ") format='%lf'" << " dx = " << in.DetecDTheta << " dy = " 
                 << in.DetecDPhi << " origin=(" << in.DetecThetaStart << "," 
                 << in.DetecPhiStart << ") w image";
    opl << "\n -------------------------------------------------------- \n";   
    opl.close();
 }
}

void LEADS_Dynamics ::WriteLaserPulse()
{
  rowvec E,B;
  rowvec r(3,fill::zeros);
  mat laser_profile(0,4); 
  
  int rt = 0;
  for(double tr = -2*in.Tau_FWHM; tr< 2*in.Tau_FWHM; tr+=dTau*in.WriteTrajecAfter)
  {
   laser_profile.resize(rt+1,4);
   la.LaserProfile(tr,r,E,B); 
   
   laser_profile(rt,0) = tr; 
   laser_profile(rt,1) = E(0); 
   laser_profile(rt,2) = E(1); 
   laser_profile(rt,3) = E(2); 
    
   rt++;
  }
 
 Print("Saving temporal profile of laser pulse"); 
 ofstream laser("./Data/laser_pulse.dat",ios::out);
 laser << "#t,Ex,Ey,Ez" << endl;
 laser_profile.save(laser,raw_ascii);
 laser_profile.clear();
 
 WritePlotInfo("laser",0,0);
 
 if(in.WriteSpatialFieldsProfile == true)
 {
  if(in.xy + in.xz + in.yz > 1 || in.xy + in.xz + in.yz < 1 )
  {Print("Select only single field snapshot plane");}
  else
  {
   if(in.xy > 0){ Print("Storing fields in xy plane");}
   if(in.xz > 0){ Print("Storing fields in xz plane");} 
   if(in.yz > 0){ Print("Storing fields in yz plane");} 
   WriteSpatialFields();
  } 
  
  if(in.RunFullSimulation == false)
  {
   VPrint("RunFullSimulation","false");
   Print("Not running full simulation");
   cout << endl << "========================================================\n";  
   exit(0);
  } 
 }
} 

void LEADS_Dynamics :: WriteSpatialFields()
{
 VPrint("Storing spatial field profiles at t",in.at_time); 
 system("mkdir -p ./Data/fields");

 double k  = 2 *pi /in.Wavelength;
 double w0 = in.BeamWaist * um;
 
 rowvec E,B;
 int N = 200;
 mat f_ex(N,N);
 mat f_ey(N,N);
 mat f_ez(N,N);
 mat f_bx(N,N);
 mat f_by(N,N);
 mat f_bz(N,N);
   
 vec x = k*w0* linspace(-4,4,N);
 vec y = k*w0* linspace(-4,4,N);
 vec z = k*w0* linspace(-40,40,N); 
  
 for(int i = 0;i< N ;i++)
 for(int j = 0;j< N ;j++)
 {
  rowvec r(3,fill::zeros); 
  if(in.xy > 0){r(0) = x(i);r(1) = y(j);r(2) = 0;}  
  if(in.xz > 0){r(0) = x(i);r(1) = 0;r(2) = z(j);}  
  if(in.yz > 0){r(0) = 0;r(1) = y(i);r(2) = z(j);}  
 
  la.LaserProfile(in.at_time,r,E,B);
  
  f_ex(i,j) = E(0);
  f_ey(i,j) = E(1);
  f_ez(i,j) = E(2);
  f_bx(i,j) = B(0);
  f_by(i,j) = B(1);
  f_bz(i,j) = B(2); 
 }
  
 f_ex.save("./Data/fields/ex.dat",raw_ascii); 
 f_ey.save("./Data/fields/ey.dat",raw_ascii); 
 f_ez.save("./Data/fields/ez.dat",raw_ascii); 
 f_bx.save("./Data/fields/bx.dat",raw_ascii); 
 f_by.save("./Data/fields/by.dat",raw_ascii); 
 f_bz.save("./Data/fields/bz.dat",raw_ascii); 
 
 WritePlotInfo("sfields",0,0);
}
 
void LEADS_Dynamics ::WriteSelectedData()
{
  //------------------------------------------------------------------------------------- 
  Print("Saving 0th elec. trajectory");
  mat store = vRP( span(0), span::all, span::all );
  mat stGama0 = vBETA( span(0), span::all, span::all );
  stGama0 = stGama0.t();
  store = store.t();  
  store.insert_cols(0,TR);
  store.insert_cols(4,stGama0);
  vec gama1(store.n_rows,fill::zeros); 
  gama1 = 1/sqrt(1 - square(store.col(4)) - square(store.col(5)) - square(store.col(6)));
  store.insert_cols(7,gama1); 
  store.save("./Data/eTraj_n0.dat",raw_ascii);  
  store.clear(); stGama0.clear(); gama1.clear();
  
  WritePlotInfo("etraj",0,0);
 
 if(in.WritePhaseSpace == true)
 WritePhaseSpace();

 if(in.WriteTrajectories == true)
 WriteBunchTrajectories(); 
}

void LEADS_Dynamics ::WritePhaseSpace()
{
 Print("Storing phase-space information");
 system("mkdir -p ./Data/phasespace");
 
 int NCols;
 
 for(int i = 0;i < nT;i += in.PhaseStoreInterval)
 {
   ofstream outps(DumpPhaseSpace(int(TR(0) + i * dTau)).c_str(),ios::binary);
   mat PhaseSpace = vRP.slice(i);
   PhaseSpace.insert_cols(3,vBETA.slice(i));
   double phaseData[PhaseSpace.n_rows][PhaseSpace.n_cols];
   for(int l =0;l< PhaseSpace.n_rows;l++)
   for(int j = 0;j< PhaseSpace.n_cols;j++)
   phaseData[l][j] = PhaseSpace(l,j);
   outps.write((char *)(&phaseData),sizeof(phaseData));
   NCols = PhaseSpace.n_cols;
   PhaseSpace.clear();
   outps.close();
 }
 WritePlotInfo("phase",NCols,0); 
 Print("Phase-space stored"); 
}

void LEADS_Dynamics ::WriteBunchTrajectories()
{
 Print("Storing bunch trajectories");
 system("mkdir -p ./Data/particles");
 
 int NCols;
 
 for(int np = 0 ; np < NP; np+=in.EveryParticle)
 {
   ofstream outtj(WriteTraj(np).c_str(),ios::binary);
   mat tRP   = vRP( span(np), span::all, span::all );
   mat tBETA = vBETA( span(np), span::all, span::all );
   tRP = tRP.t();tBETA = tBETA.t();  
   tRP.insert_cols(0,TR);
   tRP.insert_cols(4,tBETA);
   double trajData[tRP.n_rows][tRP.n_cols];
   for(int l =0;l< tRP.n_rows;l++)
   for(int j = 0;j< tRP.n_cols;j++)
   trajData[l][j] = tRP(l,j);
   outtj.write((char *)(&trajData),sizeof(trajData));
   NCols = tRP.n_cols;
   tRP.clear();tBETA.clear();
   outtj.close();
 }
 
 WritePlotInfo("butraj",NCols,0); 
 Print("Bunch trajectories stored"); 
}
