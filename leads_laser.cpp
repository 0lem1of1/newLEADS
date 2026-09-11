#include "leads.h"

LEADS_Laser ::LEADS_Laser()
{
 k        = 2 *pi /in.Wavelength;
 Omega    = in.Omega;
 Tau_FWHM = in.Tau_FWHM;
 a0       = in.a0;
 
 model = int(in.PulseModel);
 profile = (in.PulseProfile);
 w0 = in.BeamWaist * um;
 pola = int(in.Polarization);
 delta  = (pola > 0 ? 1.0/sqrt(2) : 1) ;
 phase = in.PulsePhase;
 zr = k*w0*w0 / 2.0;
 eps = zr > 0? w0/zr: 0; 

 DL_Const = e*e/(2*eps0*me*c*c*in.Wavelength);
 E_Const = me*Omega*c/e;
 M_Const = me*Omega/e;
 L_Const = 1/k;
 
}

void LEADS_Laser :: LaserProfile (double t,rowvec r,rowvec &E,rowvec &B)
{
 E.zeros(3);
 B.zeros(3);
 rowvec A(3,fill::zeros);
 
 double x1,y1,z1,x,y,z;
 x1 = r(0);y1 = r(1);z1 = r(2);
 
 double eta,Deta_Dt,Deta_Dz; 
 eta = t - z1 + phase; Deta_Dt = 1;Deta_Dz = -1;
 
 double E0,B0; 
 E0 = B0 = a0;
 
 double g_eta, alpha, derivative, order;
 
 order = pow(2.0,profile); 
 alpha = profile > 0 ? pow(2.0,order)*log(2) / pow(Tau_FWHM,order) : 0;  //
 g_eta = exp(-alpha*pow( abs(eta),order)) ;
 derivative = alpha * order * pow( abs(eta),order-1) ;
 
 double feta,d_feta; // for chirping function if later would like to introduce  
 feta = d_feta = 0;
 

 double Ax,Ay,Az;    
 double Ex,Ey,Ez;
 double Bx,By,Bz;   
 
     
 double E1,E2,E3,B1,B2,B3;
 double S[20],C[20],EPS[20],RHO[20],XI2;
 double Psi,Psi0,PsiP,PsiR,PsiG;
 double rho,zeta,xi,nu,ra;
 double Econst,Bconst;
 double wa;
 
 switch(model)
 {
    case 0:
        Ex = E0 * delta * g_eta  * Deta_Dt 
                * ( (1 + d_feta) * sin(eta + feta) + derivative * cos(eta + feta) ) ;
        Ey = E0 * sqrt(1 - SQ(delta)) * g_eta * Deta_Dt  
                * ( derivative*sin(eta+feta) - (1 + d_feta)*cos(eta+feta) );
        Ez = 0 ;
        Bx = B0 * sqrt(1 - SQ(delta)) * g_eta * Deta_Dz 
                * ( derivative * sin(eta+feta) - (1 + d_feta)*cos(eta+feta) );
        By = -B0 * delta * g_eta * Deta_Dz 
                 * ( (1 + d_feta)*sin(eta+feta) + derivative * cos(eta+feta) ) ;
        Bz = 0;

        Ax = E0 * g_eta * delta * cos(eta + feta);
        Ay = E0 * g_eta * sqrt(1 - SQ(delta)) * sin(eta + feta);
        Az = 0;

        break;
    
    case 1:

 	if(pola == 1)
        Print("For this pulse CP is not available, using LP");
      
        x = x1 * L_Const; y = y1 * L_Const; z = z1 * L_Const;
        wa = w0 *sqrt(1 + SQ(z/zr));

        ra = sqrt(x*x + y*y) ;
        xi = x  / w0; // to SI units
        nu = y  / w0;
        zeta = z / zr;
        rho = ra / w0;
        PsiP = eta + feta ;
        PsiG = atan(zeta);
        PsiR = 0.5 * k*ra*ra/R(z);
        Psi0 = 0;
        Psi = Psi0 + PsiP - PsiR + PsiG;
        Econst = E0 * E_Const * w0 * g_eta * exp(- SQ(ra/wa)) / wa;
        Bconst = Econst/c;
        XI2 = xi*xi;

        for(int i=0; i<=15; i++)
        {
            S[i] = Sn(i,z,Psi,PsiG);
            C[i] = Cn(i,z,Psi,PsiG);
            EPS[i] = PO(eps,i*1.0);
            RHO[i] = PO(rho,i*1.0);
        }
  
        E1 = Econst* ( S[0] + EPS[2]*(XI2*S[2] - RHO[4]*S[3]/4.0) + EPS[4]*(S[2]/8.0 - RHO[2]*S[3]/4.0
                       - RHO[4]*S[4]/16.0 + RHO[2]*S[4]*XI2 - RHO[6]*S[5]/8.0 - RHO[4]*S[5]*XI2/4.0 
                       + RHO[8]*S[6]/32.0)) ;
        Ex = E1/E_Const;

        E2 = Econst * xi * nu * (EPS[2]*S[2] + EPS[4]*(RHO[2]*S[4] - RHO[4]*S[5]/4.0)) ;
        Ey = E2/E_Const;

        E3 = Econst* xi * (EPS[1]*C[1] + EPS[3]*(-C[2]/2.0 + RHO[2]*C[3] - RHO[4]*C[4]/4.0) 
                   + EPS[5]*(-3.0*C[3]/8.0 - 3.0*RHO[2]*C[4]/8.0 + 17.0*RHO[4]*C[5]/16.0 
                   - 3.0*RHO[6]*C[6]/8.0 + RHO[8]*C[7]/32.0)) ;
        Ez = E3/E_Const;

        B1 = 0.0;
        Bx = B1/M_Const;

        B2 = Bconst * (S[0] + EPS[2]*(RHO[2]*S[2]/2.0 - RHO[4]*S[3]/4.0) + EPS[4]*(-S[2]/8.0 + RHO[2]*S[3]/4.0
                      + 5.0*RHO[4]*S[4]/16.0 - RHO[6]*S[5]/4.0 + RHO[8]*S[6]/32.0));
        By = B2/M_Const;
        
        B3 = Bconst * nu * (EPS[1]*C[1] + EPS[3]*(C[2]/2.0 + RHO[2]*C[3]/2.0 - RHO[4]*C[4]/4.0) 
                    + EPS[5]*(3.0*C[3]/8.0 + 3*RHO[2]*C[4]/8.0 + 3.0*RHO[4]*C[5]/16.0 - RHO[6]*C[6]/4.0 
                    + RHO[8]*C[7]/32.0));
        Bz = B3/M_Const;
  
        Ax = 0;
        Ay = 0;
        Az = 0;

        break ;
        
    case 2:
        {
        x = x1; y = y1; z = z1;
        
        double z0 = k*zr;
        double t0 = phase;
        double zeta2 = pola;
        double T = Tau_FWHM/sqrt(8*log(2));
        Complex Rc = sqrt(x*x + y*y + (z + iota0*z0)*(z + iota0*z0));
        Rc = imag(Rc) < 0 ? -Rc : Rc; 
        Complex tc = (t - t0 + iota0*z0);
        Complex tau = (tc - Rc);
        Complex p0 = (z0 * E0)/sqrt((1 - 1/z0 + 1/(T*T) + 1/(z0*z0)) * (1 - 1/z0 + 1/(T*T)));
        double phi0 = 0;	
        Complex phase0 = exp(iota0*(tau + phi0)) ;
        Complex P_0 = p0*exp(-0.5*(tau*tau/T/T)) * phase0;
        Complex Zc = z + iota0*z0;

        Complex f = (1.0 + (iota0*tau)/(T*T))*(1.0 + (iota0*tau)/(T*T)) 
                    - (1.0/(Rc*Rc))*(1.0 - tc*Rc/(T*T) + iota0*Rc);
        Complex g = -f + (2.0/(Rc*Rc))*(1.0 - tau*Rc/(T*T) + iota0*Rc);
        Complex hf = f + 1.0/(Rc*Rc);

        Complex dummyEx = (P_0)/(Rc*sqrt(1.0+zeta2*zeta2))*(f + g*x*(x + iota0*zeta2*y)/(Rc*Rc));
        Complex dummyEy = (P_0)/(Rc*sqrt(1.0+zeta2*zeta2))*(iota0*f*zeta2 + g*y*(x + iota0*zeta2*y)/(Rc*Rc));
        Complex dummyEz = (P_0)/(Rc*sqrt(1.0+zeta2*zeta2))*(g*Zc*(x + iota0*zeta2*y)/(Rc*Rc));

        double factor = 1; //* 4*PI*1e-7; initially the factor was multiplied by u0 but to convert H into B, u0 was removed

        Complex dummyBx = (P_0*hf)/(Rc*Rc*sqrt(1.0+zeta2*zeta2)*factor)*(-iota0*zeta2*Zc);
        Complex dummyBy = (P_0*hf)/(Rc*Rc*sqrt(1.0+zeta2*zeta2)*factor)*(Zc);
        Complex dummyBz = (P_0*hf)/(Rc*Rc*sqrt(1.0+zeta2*zeta2)*factor)*(iota0*zeta2*x - y);

        Ex = real(dummyEx);
        Ey = real(dummyEy);
        Ez = real(dummyEz);
        Bx = real(dummyBx);
        By = real(dummyBy);
        Bz = real(dummyBz);
        Ax = 0;
        Ay = 0;
        Az = 0;
        }
        break;        
 
   case 3:

 	if(pola == 1)
        Print("For this pulse CP is not available, using LP");
        
        x = x1 * L_Const; y = y1 * L_Const; z = z1 * L_Const;
        wa = w0 *sqrt(1 + SQ(z/zr));

        ra = sqrt(x*x + y*y) ;
        xi = x  / w0; // to SI units
        nu = y  / w0;
        zeta = z / zr;
        rho = ra / w0;
        PsiP = eta + feta ;
        PsiG = atan(zeta);
        PsiR = 0.5 * k*ra*ra/R(z);
        Psi0 = 0;
        Psi = Psi0 + PsiP - PsiR + PsiG;
        Econst = E0 * E_Const * g_eta * exp(- SQ(ra/wa)) ;
        Bconst = Econst/c;
        XI2 = xi*xi;

        for(int i=0; i<=15; i++)
        {
            S[i] = Sn(i,z,Psi,PsiG);
            C[i] = Cn(i,z,Psi,PsiG);
            EPS[i] = PO(eps,i*1.0);
            RHO[i] = PO(rho,i*1.0);
        }
        
         E1 = Econst * (  EPS[1]*RHO[1]*C[2] 
                        + EPS[3]* (-0.5*RHO[1]*C[3] + RHO[3]*C[4]-0.25*RHO[5]*C[5]) 
                        + EPS[5]* (-0.375*RHO[1]*C[4]-0.375*RHO[3]*C[5]+1.0625*RHO[5]*C[6] 
                                   -0.375*RHO[7]*C[7]+0.03125*RHO[9]*C[8])); // Er 
         Ex = (E1/E_Const) * cos(atan2(nu,xi)); 

         E2 = E1;
         Ey = (E2/E_Const) * sin(atan2(nu,xi)); 

         E3 = Econst * ( EPS[2]*(S[2]-RHO[2]*S[3])    
                       + EPS[4]*(0.5*S[3]+0.5*RHO[2]*S[4]-1.25*RHO[4]*S[5]+0.25*RHO[6]*S[6])); // Ez
         Ez = E3/E_Const;

         B1 = Bconst * ( EPS[1]*RHO[1]*C[2] 
                       + EPS[3]* (0.5*RHO[1]*C[3] + 0.5*RHO[3]*C[4]-0.25*RHO[5]*C[5])
                       + EPS[5]* (0.375*RHO[1]*C[4]+0.375*RHO[3]*C[5]+0.1875*RHO[5]*C[6] 
                                   -0.25*RHO[7]*C[7]+0.03125*RHO[9]*C[8])); // Btheta
         Bx = -(B1/M_Const) * sin(atan2(nu,xi));

         B2 =  B1; 
         By = (B2/M_Const) * cos(atan2(nu,xi));;
        
         B3 = 0; // Bz
         Bz = B3/M_Const;

        Ax = 0;
        Ay = 0;
        Az = 0;

        break ;
 }
 
 E(0) = Ex; E(1) = Ey; E(2) = Ez;
 B(0) = Bx; B(1) = By; B(2) = Bz;
 A(0) = Ax; A(1) = Ay; A(2) = Az;

}

void LEADS_Laser :: FillParams(LeadsParams &p)
{
 p.k = k; p.Omega = Omega; p.Tau_FWHM = Tau_FWHM; p.a0 = a0; p.w0 = w0;
 p.delta = delta; p.phase = phase; p.zr = zr; p.eps = eps; p.profile = profile;
 p.DL_Const = DL_Const; p.E_Const = E_Const; p.M_Const = M_Const; p.L_Const = L_Const;

 // Same expressions LaserProfile() uses inline, so the sidecar always
 // matches what the field calculation actually used.
 p.order = pow(2.0,profile);
 p.alpha = profile > 0 ? pow(2.0,p.order)*log(2) / pow(Tau_FWHM,p.order) : 0;

 p.model = model;
 p.pola = pola;
}
