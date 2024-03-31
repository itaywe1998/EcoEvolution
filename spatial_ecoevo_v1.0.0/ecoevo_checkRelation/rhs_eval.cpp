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

// [[Rcpp::export]]
double smoothstep(double x) {
  double y;
  if (x<0.0) y=0.0;
  else if (x>1.0) y=1.0;
  else y=x*x*x*(10.0+x*(-15.0+6.0*x));
  return(y);
}

// [[Rcpp::export]]
double Temp(double t,double Tmin,double C, double tE) {
  return Tmin+C*smoothstep(t/tE);
}

// [[Rcpp::export]]
List eqs(double time, NumericVector state, List pars) {
  // Parameters
  double nmin=pars["nmin"];
  double aw=pars["aw"], bw=pars["bw"], kappa=pars["kappa"];
  double Tmin=pars["Tmin"],tE=pars["tE"], C=pars["C"];
  double V=pars["s"],  rho=pars["rho"];
  double w, sw, ef, b, g, q, h2,T;
  double n, m;
  NumericVector dvdt(2);
  
  n = state[0];
  m = state[1];
  
  // if( n <nmin ){
  //   cout<<"Population died"<<endl;
  //   cout<<n<<endl;
  //   throw range_error(to_string(time));
  // }
  T=Temp(time,  Tmin, C, tE); // Vector of temperatures
  w=bw-aw*m;
  sw=w*w+V;
  ef=rho*exp(-(T-m)*(T-m)/(2.0*sw))/sqrt(sw);
  b=ef-kappa;
  g=ef*V*(T-m)/sw;
  q=V*smoothstep(n/nmin);
  h2=q/(q+V);

  dvdt[0]=(n*b)*smoothstep(n/nmin);
  dvdt[1]=h2*g;
  return(List::create(dvdt));
}
