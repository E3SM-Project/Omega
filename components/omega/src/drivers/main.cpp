
//////////////////////////////////////////////////////////////////////////////////////////
// miniWeather
// Author: Matt Norman <normanmr@ornl.gov>  , Oak Ridge National Laboratory
// This code simulates dry, stratified, compressible, non-hydrostatic fluid flows
// For documentation, please see the attached documentation in the "documentation" folder
//
//////////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <iostream>
#include <math.h>

#include "logging.h"

double _main(int argc, char **argv, double * mass, double * te);

///////////////////////////////////////////////////////////////////////////////////////
// THE MAIN PROGRAM STARTS HERE
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

  OMEGA_LOG_INFO("Starting Omega...");

  int retval = 1;
 
  double mass = -1.0;
  double te = -1.0;

  _main(argc, argv, &mass, &te);

  if (isnan(mass)) {
    std::cout << "Mass change is NaN" << std::endl; 

  } else if (abs(mass) >= 1e-9) {
    std::cout << "ERROR: Mass change magnitude is too large" << std::endl; 
    std::cout << "Change in mass: " << mass << std::endl; 
    std::cout << "Mass toleranc: " << 1e-9 << std::endl; 
   
  } else {
  
    if (isnan(te)) {
      std::cout << "Total energy change is NaN" << std::endl; 
  
    } else if (te >= 0) {
      std::cout << "ERROR: Total energy change must be negative" << std::endl; 
   
    } else if (abs(te) >= 4.5e-5) {
      std::cout << "ERROR: Total energy change magnitude is too large" << std::endl; 
      std::cout << "Change in total energy: " << te << std::endl; 
      std::cout << "Total energy toleranc: " << 4.5e-5 << std::endl; 
     
    } else {
      retval = 0;
    }
  }

  OMEGA_LOG_INFO("Omega is finished.");

  return retval;

}
