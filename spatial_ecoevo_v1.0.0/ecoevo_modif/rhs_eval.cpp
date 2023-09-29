/*
 Copyright (C) 2021 György Barabás
 This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
 you are welcome to redistribute it under certain conditions. for details,
 see the GNU General Public License Agreement (in the file COPYING.txt).
 */
#include <Rcpp.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <cmath>
#include<sstream>
#include<vector>
using namespace Rcpp;
using namespace std;
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

// [[Rcpp::export]]
NumericVector string2Vec(string str) {
  int num;
  NumericVector vec;
  stringstream ss(str);
  while ( ss >> num ) { 
    vec.push_back(num);
  }   
  return vec;
}

// [[Rcpp::export]]
string vec2String(NumericVector v) {
  string str;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing spac
  return s;
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
NumericVector Temp(NumericVector x, double t, double tE,
                   double Cmax, double Cmin, double Tmax, double Tmin) {
  return((Tmax-Tmin)*x+Tmin+((Cmin-Cmax)*x+Cmax)*smoothstep(t/tE));
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
 - state: Vector of state variables, with 2*MAX_S*L entries, where MAX_S is the   
 maximum number of species and L the number of patches. The first S*L entries are the
 densities, the second S*L entries are the trait means.
 - pars: Model parameters, given as members of a list
 Output:
 - The derivatives of the densities and trait means, as a vector in a list */
// [[Rcpp::export]]
List eqs(double time, NumericVector state, List pars) {
  //cout <<"time is "<<time<<endl;
  // Parameters
  int L=pars["L"],P = pars["P"],MAX_S=pars["MAX_S"],step=pars["step"];
  double split=pars["split"],speth=pars["speth"], eta=pars["eta"], nmin=pars["nmin"], venv=pars["venv"], tE=pars["tE"], Cmax=pars["Cmax"], Cmin=pars["Cmin"], Tmax=pars["Tmax"], Tmin=pars["Tmin"], aw=pars["aw"], bw=pars["bw"], kappa=pars["kappa"];
  NumericVector d=pars["d"], V=pars["V"], Th=pars["Th"], rho=pars["rho"], arate=pars["arate"], eps=pars["eps"];
  NumericMatrix vmat=pars["vmat"], W=pars["W"], mig=pars["mig"], a=pars["a"];
  String model=pars["model"];
  bool iteration=pars["iteration"];
  //External 
  int S=atoi(getenv("S"));
  NumericVector origins = string2Vec(getenv("origins"));
  // S-sized Matrices
  NumericMatrix F(S,S), alpha(S,S), beta(S,S);
  NumericMatrix n(S,L), m(S,L);
  // Variables
  int i, j, k, l, i_p, j_p, n_i, SR = S, newS = S, oldS = S, o,h;
  double sumgr, summig, w, sw, ef, b, bsumgr, bsummig, g, q, Omega, dm, h2,curr_m;
  NumericVector dvdt(2*MAX_S*L), x(L), T(L); 
  for (k=0; k<L; k++) x[k]=k/((double)L-1.0); // x initialize
  T=Temp(x, time, tE, Cmax, Cmin, Tmax, Tmin); // Vector of temperatures;
  for (i=0; i<S; i++) {
    for (k=0; k<L; k++){
      m(i,k)=state[MAX_S*L+k+i*L];
      n(i,k)=state[k+i*L];
      if (n(i,k)<1.0e-10) n(i,k)=0.0; // Extinction threshold
    }
  } // m & n initialize
  // Assign competition coeffs alpha_ij^k and selection pressures beta_ij^k
  for (k=0; k<L; k++) {
    // If we have temperature-dependent competition:
    if ((model=="Tdep") || (model=="Tdep_trophic")) {
      for (i=0; i<(SR-1); i++) {
        i_p = i % P;
        alpha(i,i)=eta/sqrt(2.0*V[i_p]+2.0*V[i_p]+eta*eta);
        for (j=i+1; j<SR; j++) {
          j_p = j%P;
          Omega=2.0*V[i_p]+2.0*V[j_p]+eta*eta;
          dm=m(j,k)-m(i,k);
          alpha(i,j)=eta*exp(-dm*dm/Omega)/sqrt(Omega);
          alpha(j,i)=alpha(i,j);
          beta(i,j)=2.0*V[i_p]*alpha(i,j)*dm/Omega;
          beta(j,i)=-beta(i,j)*V[j_p]/V[i_p];
        }
      }
      alpha(SR-1,SR-1)=eta/sqrt(2.0*V[SR-1]+2.0*V[SR-1]+eta*eta);
    } else { // If no temperature-dependent competition, it's much simpler:
      alpha=a;
    }
    F=funcresp(n(_,k), Th, arate, W); // Feeding rate of species i on j in patch k
    for (i=0; i<S; i++) {
      i_p = i % P;
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
        bsummig+=mig(k,l)*n(i,l)*(m(i,l)-m(i,k))/(n(i,k)+1.0e-10);
      }
      // Growth terms in the equations. V-genetic, d -dispersal
      summig*=d[i_p];
      bsummig*=d[i_p];
      w=bw-aw*m(i,k);
      sw=w*w+V[i_p];
      ef=rho[i_p]*exp(-(T[k]-m(i,k))*(T[k]-m(i,k))/(2.0*sw))/sqrt(sw);
      b=ef-kappa;
      g=ef*V[i_p]*(T[k]-m(i,k))/sw;
      q=vmat(i_p,k)*smoothstep(n(i,k)/nmin);
      h2=q/(q+venv); // Heritability
      // Assign calculated rates to vector of derivatives for output
      dvdt[k+i*L]=(n(i,k)*b+sumgr)*smoothstep(n(i,k)/1.0e-6)+summig;// dn/dt
      dvdt[MAX_S*L+k+i*L]=h2*(g-bsumgr+bsummig);// dm/dt
      // Periodic boundary conditions
      if (k==0) {
        dvdt[i*L]+=d[i]*(mig(0,1)*n(i,1)-mig(1,0)*n(i,0));
        dvdt[MAX_S*L+i*L]+=d[i]*h2*mig(0,1)*n(i,1)*(m(i,1)-m(i,0))/(n(i,0)+1.0e-10);
      }
      if (k==(L-1)) {
        dvdt[k+i*L]+=d[i]*(mig(k,k-1)*n(i,k-1)-mig(k-1,k)*n(i,k)); 
        dvdt[MAX_S*L+k+i*L]+=d[i]*h2*mig(k,k-1)*n(i,k-1)*
          (m(i,k-1)-m(i,k))/(n(i,k)+1.0e-10); 
      }
    }
  }
  //Determine the new S (Speciation Check)
  h = 4; // harsh the Speciation Condition
  for (i=0; i<S; i++) {
    o = origins[i]-1;
    for (k=0; k<L; k++){
      curr_m=m(i,k);
      if (h*abs(curr_m-T[k]) < abs(curr_m-T[o]) && n(i,k)>speth){ //if closer in trait to local env
        newS++;
      }
    }
  }
  if(newS>=MAX_S){
    throw("Surpassed Max!");
  } //Throw Exception!
  S = newS;
  if(S > oldS){
    cout<<"S is "<<S<<" oldS is "<<oldS<<endl;
  //  cout<<setprecision(4)<<n<<endl;
  //  cout<<setprecision(4)<<m<<endl;
    n_i = oldS;
    for (i=0; i<oldS; i++) {
      o = origins[i]-1;
      for (k=0; k<L; k++){
        curr_m=m(i,k);
        if (h*abs(curr_m-T[k]) < abs(curr_m-T[o]) && n(i,k)>speth){
          if(n_i>S){
            throw("Oh my god, too many aliens!");
          } //Exception Thrown
          origins.push_back(k+1);
          cout << "Species #"<<i+1<<" gave birth to the new species #"<<
            n_i+1<<" because "<<curr_m<<" is closer to the local temp "<<T[k]<<
              " of cell #"<<k+1<<" than the origin cell #"<<o+1<<" with temp "<<T[o]<<".Time:"<<time<<".Density:"<<n(i,k)<<endl;
          //Let's hope this will work, wanted behavior in comments:
          dvdt[k+n_i*L]=n(i,k)*split / step;  //newN : 0 -> oldN*split
          for(int k=0;k<L;k++){
            dvdt[MAX_S*L+k+n_i*L]=curr_m /step; //newM : 0 ->oldM - needs to be written to all cells! Most of them has n(n_i,k)=0!
          }
          dvdt[k+i*L]=-n(i,k)*split / step;   //oldN : oldN ->oldN*(1-split)
          dvdt[MAX_S*L+k+i*L]=-(curr_m-m(i,o))/step; //oldM : oldM -> M in old origin
          //I am doing all in my powers to not re-do an Speciation Event due to labeling technicalities
          n_i++;  
        }
      }
    }
    setenv("S",to_string(S).c_str(),1);
    setenv("origins",vec2String(origins).c_str(),1);
  }
  if(iteration) {
    return(List::create((state+dvdt)/step));
  }
  else return(List::create(dvdt));
}
