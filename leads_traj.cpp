#include "leads.h"

void LEADS_Dynamics ::CalForceMovePart(double tr,int i_p)
{
  /*
   Boris Algorithm : Computing The Particle Paths In An Open-Trap Sharp-Point Geometry,
   D. S. Filippychev, Computational Mathematics and Modeling,vol. 12,193–210 (2001).
   https://link.springer.com/article/10.1023/A:1012589205469
  */
  rowvec rp,beta,pbeta,bdot,B_old,E_old;
  rp = posi.row(i_p);beta = velo.row(i_p);pbeta = momt.row(i_p);bdot = accl.row(i_p);
  B_old = bfield.row(i_p);E_old = efield.row(i_p);
  
  rowvec E,B;
  la.LaserProfile(tr,rp,E,B);
  
  rowvec P,Pold,Vold;
  P = pbeta;Pold = P;Vold = beta; 
  
  double gamma = sqrt(1 + cdot(P,P));
  double gamaOld = gamma;
  double del_s = (q_part/m_part) * dTau/2.0;
  double del = del_s/gamma;
  double del2 = del*del;
  double B2 = cdot(B,B);
  double term1 = (1 - del2*B2) / (1 + del2*B2);
  double term2 =  2*del / (1 + del2*B2);
  rowvec UM = P + del_s * E ;
  rowvec UxB = Cross(UM,B);
  double udb = cdot(UM,B);
  rowvec UP = term1 * UM + term2 * UxB + term2 * del * udb * B;
  P = UP + del_s * E;
 
 if(RRFLAG != 0 && in.RRModel == 0 )
 {
    /* 
     Landau-Lifshitz Model, without derivative term
     Radiation reaction effects on radiation pressure acceleration
     M Tamburini et. al., New J. Phys. 12, 123005,(2010)
     https://doi.org/10.1088/1367-2630/12/12/123005
    */
    rowvec V = Vold; 
    rowvec FL = q_part * (E + Cross(V,B) );     
    double FL2 = cdot(FL,FL);
    double VdotE = cdot(V,E);
    rowvec FR = -RR_CONST_LL * (Cross(FL,B) + q_part * VdotE * E + SQ(gamaOld) * (FL2 - SQ(VdotE)) * V ); 
       
    pbeta = P + dTau * FR;
    gamma = sqrt(1+cdot(pbeta,pbeta));
    beta = pbeta / gamma ;
    rp = rp + beta * dTau;
    bdot = (beta - Vold)/dTau;
    
    egama(i_p) = gamma;
    posi.row(i_p) = rp; velo.row(i_p)=beta;momt.row(i_p)=pbeta;accl.row(i_p)=bdot;
    efield.row(i_p) = E;bfield.row(i_p) = B; 
 } 
 else if(RRFLAG != 0 && in.RRModel == 1 )
 {
    /* 
     Landau-Lifshitz Model, with EM field time derivative terms
     Radiation reaction effects on radiation pressure acceleration
     M Tamburini et. al., New J. Phys. 12, 123005,(2010)
     https://doi.org/10.1088/1367-2630/12/12/123005
    */
    rowvec V = Vold; 
    rowvec FL = q_part * (E + Cross(V,B) );     
    double FL2 = cdot(FL,FL);
    double VdotE = cdot(V,E);
    rowvec FR = -RR_CONST_LL * (Cross(FL,B) + q_part * VdotE * E + SQ(gamaOld) * (FL2 - SQ(VdotE)) * V ); 
    
    /*
    double gama_temp = sqrt(1+ cdot(P,P));
    rowvec P1 = P; 
    rowvec V1 = P/gama_temp;
    rowvec dR = {dTau*(V1(0)-V(0)),dTau*(V1(1)-V(1)),dTau*(V1(2)-V(2))}; 
    dR.clean(1e-6);
    rowvec dE = E - E_old;
    rowvec dB = B - B_old;
    rowvec v_grad_E = {V(0)*dE(0)/dR(0)+V(0)*dE(1)/dR(0)+V(0)*dE(2)/dR(0),
                       V(1)*dE(0)/dR(1)+V(1)*dE(1)/dR(1)+V(1)*dE(2)/dR(1),
                       V(2)*dE(0)/dR(2)+V(2)*dE(1)/dR(2)+V(2)*dE(2)/dR(2)};  
    rowvec v_grad_B = {V(0)*dB(0)/dR(0)+V(0)*dB(1)/dR(0)+V(0)*dB(2)/dR(0),
                       V(1)*dB(0)/dR(1)+V(1)*dB(1)/dR(1)+V(1)*dB(2)/dR(1),
                       V(2)*dB(0)/dR(2)+V(2)*dB(1)/dR(2)+V(2)*dB(2)/dR(2)};  
    */
     
    rowvec dE_dt = (E - E_old)/dTau; // + v_grad_E ; neglected
    rowvec dB_dt = (B - B_old)/dTau; // + v_grad_B ; neglected
    rowvec FRTotal = FR + q_part * gamaOld * RR_CONST_LL * (dE_dt + Cross(V,dB_dt));
       
    pbeta = P + dTau * FRTotal;
    gamma = sqrt(1+cdot(pbeta,pbeta));
    beta = pbeta / gamma ;
    rp = rp + beta * dTau;
    bdot = (beta - Vold)/dTau;
    
    egama(i_p) = gamma;
    posi.row(i_p) = rp; velo.row(i_p)=beta;momt.row(i_p)=pbeta;accl.row(i_p)=bdot;
    efield.row(i_p) = E;bfield.row(i_p) = B;
 } 
 else if(RRFLAG != 0 && in.RRModel == 2 )
 {
    /*
     Dynamics of emitting electrons in strong laser fields
     I V Sokolov et. al., Phys. Plasmas 16, 093115, (2009)
     https://doi.org/10.1063/1.3236748
    */
    rowvec V = Vold; 
    rowvec FL = q_part * (E + Cross(V,B) );     
    double dU_term0 = RR_CONST_SOKOLOV/(m_part + RR_CONST_SOKOLOV * cdot(V,FL)) ;
    rowvec dU = dU_term0 * (FL - cdot(V,FL)*V);	
    rowvec FR = q_part * Cross(dU,B) - SQ(gamaOld) * cdot(dU,FL) * V; 
    
    pbeta = P + dTau * FR;
    gamma = sqrt(1+cdot(pbeta,pbeta));
    beta = pbeta / gamma + dU ;
    rp = rp + beta * dTau;
    bdot = (beta - Vold)/dTau;
    
    egama(i_p) = gamma;
    posi.row(i_p) = rp; velo.row(i_p)=beta;momt.row(i_p)=pbeta;accl.row(i_p)=bdot;
    efield.row(i_p) = E;bfield.row(i_p) = B;
 }
 else if(RRFLAG != 0 && in.RRModel == 3 )
 {
    /*
     Relativistic form of Radiation Reaction
     G W Ford and R F O'Connell, Phys. Letters A 174, 182-184, (1993)
     https://doi.org/10.1016/0375-9601(93)90755-O
    */
    rowvec V = Vold; 
    rowvec FL = q_part * (E + Cross(V,B) );     
    rowvec FL_Old = q_part * (E_old + Cross(V,B_old));
    rowvec FR = RR_CONST_FORD * (gamaOld*(FL-FL_Old)/dTau - CUBE(gamaOld)* Cross(bdot,Cross(V,FL))); 
    
    pbeta = P + dTau * FR;
    gamma = sqrt(1+cdot(pbeta,pbeta));
    beta = pbeta / gamma ;
    rp = rp + beta * dTau;
    bdot = (beta - Vold)/dTau;
    
    egama(i_p) = gamma;
    posi.row(i_p) = rp; velo.row(i_p)=beta;momt.row(i_p)=pbeta;accl.row(i_p)=bdot;
    efield.row(i_p) = E;bfield.row(i_p) = B;
 }
 else
 {
  pbeta = P ;
  gamma = sqrt(1+cdot(pbeta,pbeta));
  beta = pbeta / gamma ;
  rp = rp + beta * dTau;
  bdot = (beta - Vold)/dTau;
  
  egama(i_p) = gamma;
  posi.row(i_p) = rp; velo.row(i_p)=beta;momt.row(i_p)=pbeta;accl.row(i_p)=bdot;
  efield.row(i_p) = E;bfield.row(i_p) = B;
 }
 
}


