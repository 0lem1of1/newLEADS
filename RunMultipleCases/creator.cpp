#include <iostream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
using namespace std;

int main()
{
  string oldstr = "var_value = 0.0/";
  string newstr = "var_value = ";
   
  string change0 = "" ;
  change0+= oldstr; 
  
  double val_min = 3.0;
  double val_max = 30.0;
  double dval = 3;
  
  for(double n = val_min;n<=val_max;n+= dval)
  {
   ostringstream temp; temp << n;
   string name = temp.str();
  
   system("cp -r ./NewCase junk");
   string change0 = "sed -e  \"s/" ;
   change0+= oldstr; 
   
   change0 += newstr + name 
              + "/g\" < ./junk/input.cpp > input.cpp; mv input.cpp ./junk ; mv ./junk ./val-" 
              + name; 
    
   system(change0.c_str());
  }
  
  ofstream outf("run_script",ios::out);
  for(double n = val_min;n<=val_max;n+= dval)
  {
    ostringstream temp;
    string name;
    temp << n ;
    name = temp.str();
    outf << endl << "cd val-" << name ;
    outf << endl << "make;./leads 1;cp ./Data/phasespace/time-n199 . ; ";
    outf << endl << "cd ../" << endl ;
  }
  
  for(double n = val_min;n<=val_max;n+= dval)
  {
    ostringstream temp;
    string name;
    temp << n ;
    name = temp.str();
    outf << endl << "cd val-" << name << "; basename \"$PWD\"; chmod +x ./MaxVal.gp;./MaxVal.gp;cd ..";
  }
 
  outf.close();
   
  system ("chmod +x ./run_script;") ;
  
  return 0;
}
