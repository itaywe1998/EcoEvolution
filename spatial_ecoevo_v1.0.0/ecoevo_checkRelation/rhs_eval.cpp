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


/* Apply twice continuously differentiable smoothed step function to a number x
 Input:
 - x: Distance from pole, measured in units of the pole-to-equator distance
 Output:
 - 0 if x < 0; 10*x^3-15*x^4+6*x^5 if 0 <= x <= 1; otherwise 1 */
// [[Rcpp::export]]
double smoothstep(double x) {
  double y;
  if (x<0.0) y=0.0;
  else if (x>1.0) y=1.0;
  else y=x*x*x*(10.0+x*(-15.0+6.0*x));
  return(y);
}

/* In here the range [0,1] is divided to 'cycles' amount of periods, for which 
 * symmetric rise and fall smoothstep occurs.
 * i.e for `cycles` = 5, in range [0,0.1] y = smoothstep and in [0.1,0.2] y = 1-smoothstep,
 * with the appropriate argument normalization.
 */

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


/* Temperature as a function of space (x), time (t), and some climate parameters
 Input:
 - x: Vector of distances from pole, in units of the pole-to-equator distance
 - t: Time at which temperatures are evaluated (climate change starts at t = 0)
 - tE: Time at which climate change ends (so it lasts from t = 0 to t = tE)
 - Cmax: Projected temperature increase at North Pole
 - Cmin: Projected temperature increase at equator
 - Tmax: Initial mean temperature at equator
 - Tmin: Initial mean temperature at North Pole
 Output:
 - Vector of temperatures at each location x */
// [[Rcpp::export]]
double Temp(double t, double tE,
                   double C, double T,
                   bool periodic, int cycles, bool updown) {
  if (periodic) {
    return T+C*periodic_smoothstep((t/tE),cycles, updown);
  }
  else {
    return T+C*smoothstep(t/tE);
  }
}


/* Type II functional response
 Input:
 - n: Vector of population densities of all species in a given patch
 - Th: Vector of handling times (with dummy values for resource species)
 - arate: Vector of attack rates (with dummy values for resource species)
 - W: Adjacency matrix of trophic network; W(i,j)=1 if i eats j and 0 otherwise
 Output:
 - A matrix F(i,j), the feeding rate of consumer i on resource j
 */
// all zeroes if only consumers
// [[Rcpp::export]]
NumericMatrix funcresp(NumericVector n, NumericVector Th,
                       NumericVector arate, NumericMatrix W) {
  int i, j, S=n.size();
  double Wn;
  NumericMatrix F(S,S);
  for (i=0; i<S; i++) {
    Wn=0.0;
    for (j=0; j<S; j++) Wn+=W(i,j)*n[j];
    for (j=0; j<S; j++) F(i,j)=arate[i]*W(i,j)*n[j]/(1+arate[i]*Th[i]*Wn);
  }
  return(F);
}

/* Right-hand side of dynamical equations
 Input:
 - time: Time at which function is evaluated (explicit time-dependence)
 - state: Vector of state variables, with 2*S*L entries, where S is the number
 of species and L the number of patches. The first S*L entries are the
 densities, the second S*L entries are the trait means.
 - pars: Model parameters, given as members of a list
 Output:
 - The derivatives of the densities and trait means, as a vector in a list */
// [[Rcpp::export]]
List eqs(double time, NumericVector state, List pars) {
  // Parameters
  int cycles=pars["cycles"];
  double nmin=pars["nmin"], venv=pars["venv"];
  double tE=pars["tE"], C=pars["C"],T0=pars["T0"];
  double kappa=pars["kappa"];
  double V=pars["v"], rho=pars["rho"];
  bool periodic=pars["periodic"],updown=pars["updown"];
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
