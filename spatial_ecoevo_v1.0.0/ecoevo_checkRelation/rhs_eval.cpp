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
NumericVector Temp(double t,double Tmin,double C, double tE) {
  return(Tmin+C*smoothstep(t/tE));
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
  int S=pars["S"], SR=pars["SR"], L=pars["L"];
  double eta=pars["eta"], nmin=pars["nmin"], venv=pars["venv"];
  double aw=pars["aw"], bw=pars["bw"], kappa=pars["kappa"];
  double Tmin=pars["Tmin"],tE=pars["tE"], C=pars["C"];
  NumericVector d=pars["d"], V=pars["s"], Th=pars["Th"], rho=pars["rho"];
  NumericVector arate=pars["arate"], eps=pars["eps"];
  NumericMatrix vmat=pars["vmat"], W=pars["W"], mig=pars["mig"], a=pars["a"];
  NumericMatrix T_kozai=pars["T_kozai"];
  String model=pars["model"];
  // Variables
  int i, j, k, l, ki;
  double sumgr, summig, w, sw, ef, b, bsumgr, bsummig, g, q, Omega, dm, h2;
  NumericMatrix n(S,L), m(S,L), F(S,S), alpha(S,S), beta(S,S);
  NumericVector dvdt(2*S*L), x(L), T(L);

  // Assign state variables into matrices n and m; calculate local temperatures
  for(i = 0; i < S; i++){
    for (k=0; k<L; k++) {
      n(i,k)=state[i+k*S]; // Density of species i in patch k
      if (n(i,k)<1.0e-10) n(i,k)=0.0; // Extinction threshold
      m(i,k)=state[S*L+i+k*S]; // Trait mean of species i in patch k
      
    }
  }
  cout<<"got 108"<<endl;
  
  double maxn= max(n);
  double meann = mean(n);
  double threshold = nmin;
  
  if(maxn<threshold || meann<0){
    if(maxn< threshold) cout<<"Population died"<<endl;
    if(meann<0) cout<<"Negative Profile, Simulation Broke"<<endl;
    cout<<maxn<<endl;
    throw range_error(to_string(time));
  }
  T=Temp(time,  Tmin, C, tE); // Vector of temperatures
  // Assign competition coeffs alpha_ij^k and selection pressures beta_ij^k
  cout<<"got 119"<<endl;
  
  for (k=0; k<L; k++) {
    // If we have temperature-dependent competition:
    if ((model=="Tdep") || (model=="Tdep_trophic")) {
      for (i=0; i<(SR-1); i++) {
        alpha(i,i)=eta/sqrt(2.0*V[i]+2.0*V[i]+eta*eta);
        for (j=i+1; j<SR; j++) {
          Omega=2.0*V[i]+2.0*V[j]+eta*eta;
          dm=m(j,k)-m(i,k);
          alpha(i,j)=eta*exp(-dm*dm/Omega)/sqrt(Omega);
          alpha(j,i)=alpha(i,j);
          beta(i,j)=2.0*V[i]*alpha(i,j)*dm/Omega;
          beta(j,i)=-beta(i,j)*V[j]/V[i];
        }
      }
      alpha(SR-1,SR-1)=eta/sqrt(2.0*V[SR-1]+2.0*V[SR-1]+eta*eta);
    } else { // If no temperature-dependent competition, it's much simpler:
      alpha=a;
    }
    cout<<"got 139"<<endl;
    F=funcresp(n(_,k), Th, arate, W); // Feeding rate of species i on j in patch k
    for (i=0; i<S; i++) {
      sumgr=0.0;
      bsumgr=0.0;
      // Species interaction terms in density and then trait evolution equations
      for (j=0; j<S; j++) {
        sumgr+=-n(i,k)*alpha(i,j)*n(j,k)+eps[i]*n(i,k)*F(i,j)-n(j,k)*F(j,i);
        bsumgr+=beta(i,j)*n(j,k);
      }
      summig=0.0;
      bsummig=0.0;
      // Dispersal terms in density and then trait evolution equations
      for (l=0; l<L; l++) {
        summig+=mig(k,l)*n(i,l)-n(i,k)*mig(l,k);
        bsummig+=mig(k,l)*n(i,l)*(m(i,l)-m(i,k))/(n(i,k)+nmin);
      }
      // Growth terms in the equations
      summig*=d[i];
      bsummig*=d[i];
      w=bw-aw*m(i,k);
      sw=w*w+V[i];
      ef=rho[i]*exp(-(T[k]-m(i,k))*(T[k]-m(i,k))/(2.0*sw))/sqrt(sw);
      b=ef-kappa;
      g=ef*V[i]*(T[k]-m(i,k))/sw;
      q=vmat(i,k)*smoothstep(n(i,k)/nmin);
      h2=q/(q+venv); // Heritability
      // Assign calculated rates to vector of derivatives for output
      
      dvdt[i+k*S]=(n(i,k)*b+sumgr)*smoothstep(n(i,k)/1.0e-6)+summig;
      dvdt[S*L+i+k*S]=h2*(g-bsumgr+bsummig);
      // if(k==2 && time>5e6){
      //   double expr = g+bsumgr+bsummig;
      //   double gr = g/expr, cr=bsumgr/expr, mr=bsummig;
      //   cout<<"time is " <<time<<" ratios: g-"<<g/g<<" competition- "<<bsumgr/g<<" mig-"<<bsummig/g<<endl;
      //   
      // }
    }
  }
  return(List::create(dvdt));
}
