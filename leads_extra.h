#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <vector>
#include <armadillo>
#include <cstring>
#include <sstream>
#include <chrono> 
//#include <omp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;
using namespace arma;
using namespace std::chrono; 
typedef complex<double> Complex;

template <class T> string Int2Str( const T & i ){
stringstream sstr; sstr << i; string str = sstr.str(); return str;}

template <typename T> string To_String(T const& value){
    stringstream sstr; sstr << value;return sstr.str();
}

#define pi M_PI
#define c 299792458.0
#define me  9.10938291e-31
#define e 1.60217657e-19
#define h  6.62607e-34
#define hbar  h/(2*pi)
#define re 2.8179402894E-15
#define tau0 6.26589117061524e-24
#define eps0 8.8541878e-12
#define hartee me*pow(e,3)/pow(4*pi*eps0*hbar,2) // in eV
#define um 1e-6
#define fs 1e-15
#define cons_1  1.36627476e18
#define cons_2  5.142206e11
#define iota0   Complex(0,1) 
inline double SQ(double x){return x*x ;}
inline double CUBE(double x){return x*x*x ;}
inline double PO(double x,double y){return pow(x,y);}
inline int Sign(double x){return x < 0 ? -1: 1; }
inline rowvec Cross(rowvec v,rowvec b)
{ 
 rowvec d(3,fill::zeros); 
 d(0) = v(1)*b(2) - v(2)*b(1); 
 d(1) = v(2)*b(0) - v(0)*b(2); 
 d(2) = v(0)*b(1) - v(1)*b(0); 
 return d;
}
inline vec PowVec(double x,vec a)
{
 vec b = zeros(a.n_elem);
 for(int i = 0;i< a.n_elem;i++)
 b(i) = pow(x,a(i));
 return b;
}
inline string DumpPhaseSpace(int n)
{
    stringstream result;string file;
    file =  "./Data/phasespace/time-n";
    result << n;
    return file + result.str() ;
}
inline string WriteTraj(int n)
{
    stringstream result;string file;
    file =  "./Data/particles/part-";
    result << n;
    return file + result.str() ;
} 

#define Print(msg) {char fn[]=__FILE__;fn[strlen(fn)-4] = '\0';printf("  [%s;%d]\t%s...\n",fn,__LINE__,msg);} 
#define VPrint(msg,data){char fn[]=__FILE__;fn[strlen(fn)-4]='\0';\
                        cout<<"  ["<<fn<<":"<<__LINE__<<"]\t"<<msg<<" := "<<data<<"\n";} 

 
