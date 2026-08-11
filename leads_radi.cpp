#include "leads.h"

void LEADS_Dynamics ::CalRadiation()
{
  if(in.CalculateSpectrum == true)
  {
    Print("Spectrum for range of angles");
    
    for(double th = th_st;th <= th_en;th+=dTh)
    for(double ph = ph_st;ph <= ph_en;ph+=dPh)
    {
      RadiationModule(th,ph);
      SpectrumCalculate(th,ph);
    }
   
    WritePlotInfo("spec",0,0);
    Print("Spectrum stored for ang. range"); 
  }
  
  if(in.AngularDistribution == true)
  {
    Print("Calculating angular distribution");

    int nTH = int (1+(th_en - th_st)/dTh);
    int nPH = int (1+(ph_en - ph_st)/dPh);
    vec dTheta = linspace(th_st,th_en,nTH);
    vec dPhi   = linspace(ph_st,ph_en,nPH);
     
    AngDistri = zeros(nTH,nPH);
    
    for(int i = 0;i < nTH;i++)
    for(int j = 0;j < nPH;j++)
     {
      RadiationModule(dTheta(i),dPhi(j));
      AngDistri(i,j) = TR_Integrate(dPdOmega);
     }
         
    Print("Saving angular distribution");
    AngDistri.save("./Data/angl_distri.bin",raw_binary);
    
    WritePlotInfo("angdist",nTH,nPH);
  }
  
  if(in.SignalEmitted == true)
  {
    Print("Calculating signal emitted");
     
    for(double th = th_st;th <= th_en;th+=dTh)
    for(double ph = ph_st;ph <= ph_en;ph+=dPh)
    {
      RadiationModule(th,ph);
      CalcEmittedSignal(th,ph);
    }
    
    WritePlotInfo("signal",0,0);
     
    
    Print("Emitted signal stored"); 
  
  }
 
}

void LEADS_Dynamics :: RadiationModule(double detcTheta,double detcPhi)
{
  double r_det = 1.0E+8;
  rowvec detector = {r_det*sin(detcTheta)*cos(detcPhi),r_det*sin(detcTheta)*sin(detcPhi),r_det*cos(detcTheta)}; 
    
  double modNNU2,kappa;
  rowvec R,ND,rp,vbeta,vbdot,NNU;
  dPdOmega = zeros(nT); 
  
  for(int i = 0;i< nT ; i++)
  {
    modNNU2 = 0;
  
    for(int np = 0; np < NP; np++)
    {
      rp = {vRP(np,0,i),vRP(np,1,i),vRP(np,2,i)}; 
      vbeta = {vBETA(np,0,i),vBETA(np,1,i),vBETA(np,2,i)}; 
      vbdot = {vBETADOT(np,0,i),vBETADOT(np,1,i),vBETADOT(np,2,i)}; 
      
      R = detector - rp;  
      ND = R/norm(R);
      kappa = 1 - cdot(ND,vbeta);
      
      NNU = Cross(ND, Cross( (ND - vbeta) , vbdot)) ; 
      
      vNNU(np,0,i) = NNU(0) ;
      vNNU(np,1,i) = NNU(1) ;
      vNNU(np,2,i) = NNU(2) ;
      vT(np,i) = TR(i) - as_scalar(cdot(ND,rp));// both the expressions are valid
      //vT(np,i) = TR(i) + norm(R) - r_det; // t_ret = t - |R|/c  => t = t_ret + |R|/c
                                          // r_det is substracted to keep the actual
                                          // no. down, it wont matter as anyway later
                                          // we take | | of integrand and it will be 
                                          // just a constant and | | = 1; 
      vKAPPA(np,i) = kappa;               
      modNNU2 += SQ(norm(NNU)) / pow(kappa,5); // Eq. 14.38 : Jackson
    }
     dPdOmega(i) = modNNU2; // need to integrate w.r.to TR to calculate total energy emitted at t_ret
  }
}

void LEADS_Dynamics :: SpectrumCalculate(double thSpec,double phSpec)
{
 double s_start = in.SpectrumStart;
 double s_end = in.SpectrumEnd;
 double s_ds = in.SpectrumDO;

 int nS = int (1 + (s_end-s_start)/s_ds);
 vec sVec = linspace(s_start,s_end,nS);
 vec s = in.UseLogScale > 0 ? PowVec(10,sVec) : sVec;
 vec spec_add = zeros(nS);
 
 cx_vec integrand_x(nT); integrand_x.fill(Complex(0,0));
 cx_vec integrand_y(nT); integrand_y.fill(Complex(0,0));
 cx_vec integrand_z(nT); integrand_z.fill(Complex(0,0));
 cx_vec integral(3);integral.fill(Complex(0,0));
  
 /*
  for(int nsp = 0;nsp < nS;nsp++)
  {
    integral.fill(Complex(0,0));
    for(int np = 0; np < NP; np++)
    {
      for(int i = 0;i< nT;i++)
      {
       integrand_x(i) =  vNNU(np,0,i) * exp(iota0 * s(nsp) * vT(np,i)) / SQ(vKAPPA(np,i)); // Eq. 14.62: Jackson 
       integrand_y(i) =  vNNU(np,1,i) * exp(iota0 * s(nsp) * vT(np,i)) / SQ(vKAPPA(np,i)); // Eq. 14.62: Jackson 
       integrand_z(i) =  vNNU(np,2,i) * exp(iota0 * s(nsp) * vT(np,i)) / SQ(vKAPPA(np,i)); // Eq. 14.62: Jackson 
      }
      
      integral(0) = integral(0) + SP_Integrate(integrand_x);
      integral(1) = integral(1) + SP_Integrate(integrand_y);
      integral(2) = integral(2) + SP_Integrate(integrand_z);
    }
    spec_add(nsp) =  as_scalar(real(cdot(integral,integral))); 
  }
 */
  
  
  for(int np = 0; np < NP; np++)
  {
    for(int nsp = 0;nsp < nS;nsp++)
    {
      for(int i = 0;i< nT;i++)
      {
       integrand_x(i) =  vNNU(np,0,i) * exp(iota0 * s(nsp) * vT(np,i)) / SQ(vKAPPA(np,i)); // Eq. 14.62: Jackson 
       integrand_y(i) =  vNNU(np,1,i) * exp(iota0 * s(nsp) * vT(np,i)) / SQ(vKAPPA(np,i)); // Eq. 14.62: Jackson 
       integrand_z(i) =  vNNU(np,2,i) * exp(iota0 * s(nsp) * vT(np,i)) / SQ(vKAPPA(np,i)); // Eq. 14.62: Jackson 
      }
      
      integral(0) = SP_Integrate(integrand_x);
      integral(1) = SP_Integrate(integrand_y);
      integral(2) = SP_Integrate(integrand_z);
      
      spec_add(nsp) = spec_add(nsp) + as_scalar(real(cdot(integral,integral))); 
    }
  }
 
  
  VPrint("Saving Specturm, theta",thSpec*180/pi);
  VPrint("Saving Specturm, phi",phSpec*180/pi);
  
  ofstream ospec("./Data/spectrum.dat",ios::app); 
  mat store;
  vec th_spec(s.n_rows);th_spec.fill(thSpec*180/pi); 
  vec ph_spec(s.n_rows);ph_spec.fill(phSpec*180/pi); 
  store.insert_cols(0,th_spec);store.insert_cols(1,ph_spec);
  store.insert_cols(2,s);store.insert_cols(3,spec_add);
  store.save(ospec,raw_ascii);
  store.clear();
  ospec << "\n\n"; 
  ospec.close(); 
  
  integrand_x.clear();integrand_y.clear();integrand_z.clear();
  integral.clear();s.clear();spec_add.clear();sVec.clear();
}

void LEADS_Dynamics :: CalcEmittedSignal(double thSpec,double phSpec)
{
  double r_det = 1.0E+8;
  double R;
  
  vec t_lab_np = zeros(NP);
  vec t_lab = zeros(nT);
  E_Signal = zeros(nT,3);
  //B_Signal = zeros(nT,3);
    
  for(int i = 0;i< nT ; i++)
  { 
    t_lab_np.fill(0);
    for(int np = 0; np < NP; np++)
    {
      R = vT(np,i) - TR(i) ;  
      
      if(abs(R) > 1.0E-6)
      {      
       E_Signal(i,0) =  E_Signal(i,0) + vNNU(np,0,i)/ ( R * CUBE(vKAPPA(np,i))); // Eq. 14.14: Jackson
       E_Signal(i,1) =  E_Signal(i,1) + vNNU(np,1,i)/ ( R * CUBE(vKAPPA(np,i))); // Eq. 14.14: Jackson
       E_Signal(i,2) =  E_Signal(i,2) + vNNU(np,2,i)/ ( R * CUBE(vKAPPA(np,i))); // Eq. 14.14: Jackson        
      }
      t_lab_np(np) = vT(np,i); 
    }
    
    t_lab(i) = mean(t_lab_np);
  }
  
  VPrint("Saving signal, theta",thSpec*180/pi);
  VPrint("Saving signal, phi",phSpec*180/pi);
  
  ofstream osig("./Data/emitted_signal.dat",ios::app); 
  mat store;
  vec th_spec(nT);th_spec.fill(thSpec*180/pi); 
  vec ph_spec(nT);ph_spec.fill(phSpec*180/pi); 
  store.insert_cols(0,th_spec);
  store.insert_cols(1,ph_spec);
  store.insert_cols(2,TR);
  store.insert_cols(3,t_lab);
  store.insert_cols(4,E_Signal);
  store.insert_cols(7,dPdOmega); 
  store.save(osig,raw_ascii);
  store.clear();
  osig << "\n\n"; 
  osig.close(); 
  
}

