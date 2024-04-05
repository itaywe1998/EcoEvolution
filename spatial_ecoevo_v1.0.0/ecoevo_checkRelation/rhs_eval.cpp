/*
 Copyright (C) 2021 György Barabás
 This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
 you are welcome to redistribute it under certain conditions. for details,
 see the GNU General Public License Agreement (in the file COPYING.txt).
 */

#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <string>

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
double Temp(double t,double Tmin,double C, double tE) {
  return(Tmin+C*smoothstep(t/tE));
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
  double nmin=pars["nmin"];
  double aw=pars["aw"], bw=pars["bw"], kappa=pars["kappa"];
  double Tmin=pars["Tmin"],tE=pars["tE"], C=pars["C"];
  double rho=pars["rho"],v = pars["v"];
  // Variables
  int i, j, k, l;
  double sw,w, ef, b, g, h2;
  double n,m,T;
  NumericVector dvdt(2);
  
  n=state[0]; 
  m=state[1]; 
  
  if(n<nmin){
    cout<<"Population died"<<endl;
    cout<<n<<endl;
    throw range_error(to_string(time));
  }
  T=Temp(time,  Tmin, C, tE); // Vector of temperatures
  cout<<T<<"is T ";
  // Assign competition coeffs alpha_ij^k and selection pressures beta_ij^k
      w=bw-aw*m;
      sw=w*w+v;
      ef=rho*exp(-(T-m)*(T-m)/(2.0*sw))/sqrt(sw);
      b=ef-kappa;
      cout<<ef<<" is f"<<endl;
      g=ef*v*(T-m)/sw;
      h2=0.5; // Heritability
      // Assign calculated rates to vector of derivatives for output
      
      dvdt[0]=n*b;
      dvdt[1]=h2*g;
    }
  }
  return(List::create(dvdt));
}
