//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>

#include "eos.h"
#include "eoPNJL.h"

using namespace std;
double pnjl_ed_bounds[] = {0.0, 3.0, 5.0, 40.0}; //short, low, high
double pnjl_nB_bounds[] = {0.0, 0.2, 0.8};
double pnjl_nnB_list[]  = {201,81,81};
double pnjl_ne_list[]   = {61,100,71};
int    pnjl_iitab = 0;

// ---- auxiliary EoS class. Two instances (objects) of this class will be
// created
// Huge thanks to Yurii Karpenko :)
class EoSauxPNJL {
  double emax, nmax, emin, nmin;
  int ne, nn; //number of epsilon and nb points
  double** ptab, **Ttab, **mubtab, **mustab, **stab;

 public:
  EoSauxPNJL(char* filename, int itab);
  EoSauxPNJL(const EoSauxPNJL& ee);
  ~EoSauxPNJL();
  double get_emax() { return emax; }
  double get_nmax() { return nmax; }
  double getne() { return ne; }
  double getnn() { return nn; }
  void get(double e, double nb, double& p, double& T, double& mub, double& mus);
  double p(double e, double nb);
};

EoSauxPNJL::EoSauxPNJL(char* filename, int itab) {
   ne = pnjl_ne_list[itab];
   nn = pnjl_nnB_list[itab];
   ptab   = new double* [ne];
   Ttab   = new double* [ne];
   mubtab = new double* [ne];
   mustab = new double* [ne];
   stab   = new double* [ne];
   for (int i = 0; i < ne; i++) {
     ptab[i]   = new double[nn];
     Ttab[i]   = new double[nn];
     mubtab[i] = new double[nn];
     mustab[i] = new double[nn];
     stab[i]   = new double[nn];
   }
   double* e = new double[ne];
   double* n = new double[nn];
    
   double de, dn;
   if(itab == 0){
      de = 0.05; dn = 0.001;
   }else if(itab == 1){
      de = 0.05; dn = 0.01;
   }else if(itab == 2){
      de = 0.5;  dn = 0.01;
   }else{
     de = 0.;  dn = 0.;
     cerr << "error in itab value: " << itab << " (possible values are 0, 1 or 2)" << endl;
     exit (EXIT_FAILURE);
   }
    
    char filename1[200];
    char filename2[200];
    char filename3[200];

    sprintf(filename1, "%sPNJL_Peos_%d.dat",   filename, itab);
    sprintf(filename2, "%sPNJL_Teos_%d.dat",   filename, itab);
    sprintf(filename3, "%sPNJL_Mueos_%d.dat",  filename, itab);
    ifstream f_p(filename1);
    ifstream f_T(filename2);
    ifstream f_muB(filename3);

  if (!f_p.good()) {
    cout << "I/O error with " << filename1 << endl;
      exit(1);
  }
  if (!f_T.good()){
    cout << "I/O error with " << filename2 << endl;
      exit(1);
  }
  if (!f_muB.good()){
    cout << "I/O error with " << filename3 << endl;
      exit(1);
  }
for (int in = 0; in < nn; in++){
    for (int ie = 0; ie < ne; ie++) {
        if(itab == 2) e[ie] = 5.0 + ie * de;
        else e[ie] = ie * de;
        n[in] = in * dn;
        double skip;
        f_p   >> skip >> skip >> ptab[ie][in];       // --> p[GeV/fm3]
        f_T   >> skip >> skip >> Ttab[ie][in];       // --> T[GeV]
        f_muB >> skip >> skip >> mubtab[ie][in];     // --> mub[GeV]
       // if(itab == 1) cout << "Ttab["<<ie<<"]["<<in<<"]: "<<Ttab[ie][in] << endl;
        mustab[ie][in] = 0.0;
        mubtab[ie][in] *= 3;
}}
  emin = e[0];
  emax = e[ne - 1];
  nmin = n[0];
  nmax = n[nn - 1];


  delete[] e;
  delete[] n;
}

EoSauxPNJL::EoSauxPNJL(const EoSauxPNJL& ee){
    emax = ee.emax;
    nmax = ee.nmax;
    emin = ee.emin;
    nmin = ee.nmin;
    ne   = ee.ne;
    nn   = ee.nn;
    ptab    = new double* [ne];
    Ttab    = new double* [ne];
    mubtab  = new double* [ne];
    mustab  = new double* [ne];
    stab    = new double* [ne];
    
    for (int i = 0; i < ne; i++) {
        ptab[i]   = new double[nn];
        Ttab[i]   = new double[nn];
        mubtab[i] = new double[nn];
        mustab[i] = new double[nn];
        stab[i]   = new double[nn];
    }
    for (int i = 0; i < ne; i++) {
        for (int j = 0; j < nn; j++) {
            ptab[i][j]   = ee.ptab[i][j];
            Ttab[i][j]   = ee.Ttab[i][j];
            mubtab[i][j] = ee.mubtab[i][j];
            stab[i][j]   = ee.stab[i][j];
            mustab[i][j] = ee.mustab[i][j];
    }}
}
EoSauxPNJL::~EoSauxPNJL() {
  for (int i = 0; i < ne; i++) {
    delete[] ptab[i];
    delete[] Ttab[i];
    delete[] mubtab[i];
    delete[] mustab[i];
    delete[] stab[i];
  }
  delete[] ptab;
  delete[] Ttab;
  delete[] mubtab;
  delete[] mustab;
  delete[] stab;
}

void EoSauxPNJL::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {
  if (e < 0.) {
    T = mub = mus = p = 0.;
    return;
  }
  const double de = (emax - emin) / (ne - 1);
  const double dn = (nmax - nmin) / (nn - 1);
  int ie = (int)((e  - emin) / de);
  int in = (int)((nb - nmin) / dn);
  if (ie < 0) ie = 0;
  if (in < 0) in = 0;
  if (ie > ne - 2) ie = ne - 2;
  if (in > nn - 2) in = nn - 2;
  const double em = e  - emin - ie * de;
  const double nm = nb - nmin - in * dn;

  double we[2] = {1. - em / de, em / de};
  double wn[2] = {1. - nm / dn, nm / dn};

  T = mub = mus = p = 0.0;
  for (int je = 0; je < 2; je++){
      for (int jn = 0; jn < 2; jn++) {
        p   += we[je] * wn[jn] * ptab[ie + je][in + jn];
        T   += we[je] * wn[jn] * Ttab[ie + je][in + jn];
        mub += we[je] * wn[jn] * mubtab[ie + je][in + jn];
        mus += we[je] * wn[jn] * mustab[ie + je][in + jn];
  }}
    
  if(p < 0.0) p = 0.0;
}

double EoSauxPNJL::p(double e, double nb) {
  if (e < 0.) return 0.0;
  const double de = (emax - emin) / (ne - 1);
  const double dn = (nmax - nmin) / (nn - 1);
  int ie = (int)((e - emin) / de);
  int in = (int)((nb - nmin) / dn);
  if (ie < 0) ie = 0;
  if (in < 0) in = 0;
  if (ie > ne - 2) ie = ne - 2;
  if (in > nn - 2) in = nn - 2;
  const double em = e - emin - ie * de;
  const double nm = nb - nmin - in * dn;
  double we[2] = {1. - em / de, em / de};
  double wn[2] = {1. - nm / dn, nm / dn};

  double p = 0.0;
  if(in == 0){ // on the "left" edge of the table
    for (int je = 0; je < 2; je++)
           p += we[je]  * ptab[ie + je][0];
  }else{
    for (int je = 0; je < 2; je++)
        for (int jn = 0; jn < 2; jn++)
             p += we[je] * wn[jn] * ptab[ie + je][in + jn];
  }
   
  if (e == 0. && nb == 0. ) { p = ptab[0][0];}
  if (p < 0.0) p = 0.0;
  if(!(p >=0 || p < 0)) cout << "Wrong PRESSURE: " <<  p << endl;
  return p;
}

eoPNJL::eoPNJL(char* filename) {
    for(int i = 0; i < 3; i++){
        EoSauxPNJL object(filename, i);
        eos_t.push_back(object);
    }
}

eoPNJL::~eoPNJL() {
}

void eoPNJL::eosranges(double &_emax, double &_e0, double &_nmax,
                          double &_n0, int &_ne, int &_nn) {
 _emax = eos_t[2].get_emax();
 _nmax = eos_t[2].get_nmax();
 _e0 = 0.;  // I don't think it is a meaningul param, so return 0
 _n0 = 0.;  // I don't think it is a meaningul param, so return 0
 _ne = eos_t[2].getne(); // This does not make sense either since there are two EoS
 _nn = eos_t[2].getnn(); // tables, so I return the size of larger table.
}

void eoPNJL::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
    pnjl_iitab = -99;
    
    if((e >= 0.0 && e <= 3.0) && (nb >= 0.0 && nb <= 0.2))            pnjl_iitab = 0;
        else if ((e >= 0.0 && e <= 5.0)  && (nb >= 0.0 && nb <= 0.8)) pnjl_iitab = 1;
        else if ((e > 5.0 && e <= 40.0)  && (nb >= 0.0 && nb <= 0.8)) pnjl_iitab = 2;
    
    if (pnjl_iitab == -99){
        p = 0.2964 * e;
        T = 0.15120476935 * pow(e, 0.25);
        mub = mus = 0.0;
    }
    else eos_t[pnjl_iitab].get(e, nb, p, T, mub, mus);
    muq = 0.0;
    
}

double eoPNJL::p(double e, double nb, double nq, double ns) {

    pnjl_iitab = -99;
    
    if((e >= 0.0 && e <= 3.0) && (nb >= 0.0 && nb <= 0.2))            pnjl_iitab = 0;
        else if ((e >= 0.0 && e < 5.0) && (nb >= 0.0 && nb <= 0.8))   pnjl_iitab = 1;
        else if ((e >= 5.0 && e <= 40.0) && (nb >= 0.0 && nb <= 0.8)) pnjl_iitab = 2;
    
    if(pnjl_iitab == -99) return 0.2964 * e;
    else return eos_t[pnjl_iitab].p(e, nb);
}

void eoPNJL::eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) {
    cout << " >>>>  no original table avaliable <<<< " << endl;
}
