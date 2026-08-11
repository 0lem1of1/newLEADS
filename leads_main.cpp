#include "leads.h"

int main(int argc,char **argv)
{
  if(argc < 2)
  {
    cout << endl << "Usage : ./leads NO_OF_THREADS"  << endl << endl;
    return 0;
  }
  int Num_Threads = atoi(argv[1]);
    
  auto StartTime = high_resolution_clock::now(); 
  
  LEADS_Simulation lsimu;
  lsimu.Run_Simulation(Num_Threads);
  
  auto EndTime = high_resolution_clock::now(); 
  
  auto dsec = duration_cast<seconds>(EndTime - StartTime);   
  auto dmins = duration_cast<minutes>(EndTime - StartTime);   
  auto dhrs = duration_cast<hours>(EndTime - StartTime);   
  
  cout << endl << "========================================================\n";  
  cout << " Run Time (sec/mins/hrs): " ;
  cout << "  " << dsec.count() << setw(5) << dmins.count() << setw(5) << dhrs.count(); 
  cout << endl << "========================================================\n\n";
      
  return 0;
}
