/*
 Copyright (C) 2021 György Barabás
 This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
 you are welcome to redistribute it under certain conditions. for details,
 see the GNU General Public License Agreement (in the file COPYING.txt).
*/


#include <Rcpp.h>
#include <math.h>
#include <iostream>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double smoothstep(double x) {
  double y;
  if (x<0.0) y=0.0;
  else if (x>1.0) y=1.0;
  else y=x*x*x*(10.0+x*(-15.0+6.0*x));
  return(y);
}

// [[Rcpp::export]]
double periodic_smoothstep(double x, int cycles, bool updown) {
  double y;
  if (x<0.0) y=0.0;
  else if (x>1.0) y=0.0;
  else {
    y = sin(cycles*M_PI*x);
    if (!updown) y = abs(y);        
  }
  return y;
}


// [[Rcpp::export]]
double Temp(double t, double tE,double C, double T,bool periodic, int cycles, bool updown) {
  if (periodic) {
    return T+C*periodic_smoothstep((t/tE),cycles, updown);
  }
  else {
    return T+C*smoothstep(t/tE);
  }
}

// [[Rcpp::export]]
List eqs(double time, NumericVector state, List pars) {
  // Parameters
  double nmin=pars["nmin"];
  double tE=pars["tE"], C=pars["C"],T0=pars["T0"];
  double kappa=pars["kappa"];
  double V=pars["v"], rho=pars["rho"];
  bool periodic=pars["periodic"],updown=pars["updown"];
  int cycles=pars["cycles"];
  // Variables
  double f, g,h2;
  double n, m,T;
  NumericVector dvdt(2);
  n=state[0]; // Density of species i in patch k
  m=state[1]; // Trait mean of species i in patch k
  if(n<nmin){
    cout<<time<<endl;
    throw range_error(to_string(time));
  } 
  T=Temp(time, tE, C, T0,  periodic, cycles , updown); // Vector of temperatures
  f=rho*exp(-(T-m)*(T-m)/(2.0*V))/sqrt(V);
  g=f*(T-m);
  h2=0.5; // Heritability
  // Assign calculated rates to vector of derivatives for output
  dvdt[0]=n*(f-kappa);
  dvdt[1]=h2*g;
  return(List::create(dvdt));
}
