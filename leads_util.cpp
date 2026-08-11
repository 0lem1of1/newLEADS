#include "leads.h"

void LEADS_Simulation  ::  Greetings()
{
  system("clear");
  cout << endl << "========================================================";
  cout << endl << "  LEADS - Laser Electron interAction Dynamics Simulator ";
  cout << endl << "    Amol Holkundkar, Dept. of Physics, BITS - Pilani    ";
  cout << endl << "    e-mail: amol.holkundkar@pilani.bits-pilani.ac.in    ";
  cout << endl << "========================================================";
  cout << endl;
}

void LEADS_Simulation  :: Requirements()
{
    Greetings(); 
    
    char ANS;
    int Error = 1;
    int Output;
    int DATA = 1;

    Greetings();

    ofstream ofile("./Data/test.txt",ios::out);
    cout << endl << "Checking Output Directory...............";
    if(ofile.fail())
    {
        Output = 0;
        cout << "       [FAILED]";
    }
    else
    {
        Output = 1;
        cout << "         [OK]  ";
        cout << endl << "Checking previous data .................";
        ifstream inp2("./Data/input.cpp",ios::in);

        if(inp2.fail())
        {
            DATA = 1;
            cout << "         [Ok]  ";
        }
        else
        {
            DATA = 1;
            cout << "       [Failed]  ";
        }
        inp2.close();
    }
    system("rm -f ./Data/test.txt");

    Error = Output * DATA;
    if(Error == 0)
    {
        cout << endl << endl << "========================================================";
        cout << endl <<         "     Following errors occured, aborting program.        ";
        cout << endl <<         "========================================================";
        cout << endl;

        if(Output == 0)
            cout << endl << "# Unable to find the output folder './Data/' in current directory.";
        if(DATA == 0)
            cout << endl << "# Output directory is not empty." ;
    }

    if(Error != 1 )
    {
        cout << endl << endl;
        exit(EXIT_FAILURE);
    }

    system("cp input.cpp ./Data/");
    cout << endl << endl << "========================================================";
    cout << endl <<         "         Basic requirements met, starting LEADS         ";
    cout << endl <<         "========================================================";
    cout << endl << endl;

    ofile.close();
}


double LEADS_Dynamics :: TR_Integrate(vec func1)
{
    double from = TR(0); 
    double to = TR(TR.n_elem-1);
    double dx0 = TR(1) - TR(0);
    double dxi = dx0/in.InterpolateLevel;

        
    double* x1a = &TR(0) ;
    double* func1a = &func1(0) ; // convert vector to array

    // interpolating data
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, TR.n_elem);
    gsl_spline_init (spline, x1a, func1a, TR.n_elem);
 
    int nxi = int( (to - from)/dxi );
    vec x = linspace(from,to,nxi); 
    vec func(nxi);
 
    for (int i = 0 ; i < nxi ; i++)  
    func(i) = gsl_spline_eval (spline, x(i), acc) ;
    
    gsl_spline_free (spline); gsl_interp_accel_free(acc);
    
    return as_scalar(trapz(x,func));
    
    /*
    int n = x.size();

    double h0 = (to - from) / n;
    double sum1 = 0.0;
    double sum2 = 0.0;

    if(abs(func(0)) > 100)
    func(0) = 0;

    for(int i = 0; i < n; i++)
    sum1 += func(i);

    sum2 = sum1 - func(0);

    return (h0 / 6.0) * (func(0) + func(n-1) + 4.0 * sum1 + 2.0 * sum2) ;
    */
}

Complex LEADS_Dynamics :: SP_Integrate(cx_vec func1)
{
    vec re_func1 = real(func1);
    vec im_func1 = imag(func1);
    
    double from = TR(0); 
    double to = TR(TR.n_elem-1);
    double dx0 = TR(1) - TR(0);
    double dxi = dx0/in.InterpolateLevel;
    int nxi = int(1+(to - from)/dxi);
    vec x = linspace(from,to,nxi); 
        
    double* x1a = &TR(0) ;
     
    double* func1a_re = &re_func1(0) ; // convert vector to array
    double* func1a_im = &im_func1(0) ; // convert vector to array
        
    // interpolating data
    gsl_interp_accel *acc_re = gsl_interp_accel_alloc ();
    gsl_interp_accel *acc_im = gsl_interp_accel_alloc ();
    gsl_spline *spline_re = gsl_spline_alloc (gsl_interp_cspline, TR.n_elem);
    gsl_spline *spline_im = gsl_spline_alloc (gsl_interp_cspline, TR.n_elem);
    gsl_spline_init (spline_re, x1a, func1a_re, TR.n_elem);
    gsl_spline_init (spline_im, x1a, func1a_im, TR.n_elem);
    
    vec func_re(nxi); vec func_im(nxi);
    for (int i = 0 ; i < nxi ; i++)
    {  
     func_re(i) = gsl_spline_eval (spline_re, x(i), acc_re);
     func_im(i) = gsl_spline_eval (spline_im, x(i), acc_im);
    }
    gsl_spline_free (spline_re); gsl_interp_accel_free(acc_re);
    gsl_spline_free (spline_im); gsl_interp_accel_free(acc_im);
    
    double res_re = as_scalar(trapz(x,func_re));
    double res_im = as_scalar(trapz(x,func_im));
    Complex result0 = Complex(res_re,res_im);
    
    return result0;  
 
    /*
    int n = x.size();

    double h0 = (to - from) / n;
    double sum1r = 0.0;double sum1i = 0.0;
    double sum2r = 0.0;double sum2i = 0.0;

    if(abs(func_re(0)) > 100)
    func_re(0) = 0;
     
    if(abs(func_im(0)) > 100)
    func_im(0) = 0;
 
     
    for(int i = 0; i < n; i++)
    {sum1r += func_re(i);sum1i += func_im(i);}

    sum2r = sum1r - func_re(0);
    sum2i = sum1i - func_im(0);
    
    double res_re =  (h0 / 6.0) * (func_re(0) + func_re(n-1) + 4.0 * sum1r + 2.0 * sum2r) ;
    double res_im =  (h0 / 6.0) * (func_im(0) + func_im(n-1) + 4.0 * sum1i + 2.0 * sum2i) ;
        
    Complex result0 = (res_re,res_im);
    return result0;
    */
    
}


