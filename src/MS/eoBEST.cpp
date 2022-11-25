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
#include "eoBEST.h"

using namespace std;
double ed_bounds[] = {0.0, 0.0036, 0.015, 0.045, 0.455, 20.355, 650.};
double nB_bounds[] = {0.0, 0.0025, 0.015, 0.045, 0.25,  3.5,    12.0};
double ne_list[]   = {61,  60,  61,  122, 200, 400};
double nnB_list[]  = {501, 301, 181, 251, 351, 251};
int    iitab = 0;
double hbarC = 0.19733;
// ---- auxiliary EoS class. Two instances (objects) of this class will be
// created
// Huge thanks to Yurii Karpenko :)
class EoSauxBEST {
  double emax, nmax, emin, nmin;
  int ne, nn; //number of epsilon and nb points
  double ** ptab, **Ttab, **mubtab, **mustab, **stab;
  double **epsilontab,   **nbtab;  //original table
  double tmax, mubmax, tmin, mubmin;
  int nT, nmuB;
    
 public:
  EoSauxBEST(char* filename, int itab);
  EoSauxBEST(const EoSauxBEST& ee);
  ~EoSauxBEST();
  double get_emax() { return emax; }
  double get_nmax() { return nmax; }
  double getne() { return ne; }
  double getnn() { return nn; }
  void get(double e, double nb, double& p, double& T, double& mub, double& mus);
  double p(double e, double nb);
  //from original:
  void geteosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p);
};

EoSauxBEST::EoSauxBEST(char* filename, int itab) {
  ne = ne_list[itab];
  nn = nnB_list[itab];
  ptab = new double* [ne];
  Ttab = new double* [ne];
  mubtab = new double* [ne];
  mustab = new double* [ne];
  stab = new double* [ne];
  for (int i = 0; i < ne; i++) {
    ptab[i] = new double[nn];
    Ttab[i] = new double[nn];
    mubtab[i] = new double[nn];
    mustab[i] = new double[nn];
    stab[i] = new double[nn];
  }
  double* e = new double[ne];
  double* n = new double[nn];
  double de = (ed_bounds[itab+1]  - ed_bounds[itab]) / ne;
  double dn = (nB_bounds[itab+1]  - nB_bounds[0])    / nn;

    char filename1[200];
    char filename2[200];
    char filename3[200];

    sprintf(filename1, "%sBEST_eos_p_%d.dat",   filename, itab);
    sprintf(filename2, "%sBEST_eos_T_%d.dat",   filename, itab);
    sprintf(filename3, "%sBEST_eos_muB_%d.dat", filename, itab);
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
for (int ie = 0; ie < ne; ie++) {
    for (int in = 0; in < nn; in++){
        e[ie] = ed_bounds[itab] + ie * de;
        n[in] = nB_bounds[0]    + in * dn;
        f_p >>   ptab[ie][in];       // --> p[GeV/fm3]
        f_T >>   Ttab[ie][in];       // --> T[GeV]
        f_muB >> mubtab[ie][in];     // --> mub[GeV]
        mustab[ie][in] = 0.0;
}}
  emin = e[0];
  emax = e[ne - 1];
  nmin = n[0];
  nmax = n[nn - 1];

  delete[] e;
  delete[] n;
  

// -----   READING ORIGINAL  ------
//  char parameters[200] = "Files_PAR_143_350_3_93_143_286";
  char filename4[200];
  char filename5[200];
  sprintf(filename4, "%s/BEST_EnerDens_Final.dat", filename);
  sprintf(filename5, "%s/BEST_BarDens_Final.dat",  filename);
    
  //  sprintf(filename4, "%s%s/BEST_EnerDens_Final.dat", filename,  parameters);
  //  sprintf(filename5, "%s%s/BEST_BarDens_Final.dat",  filename,  parameters);
    
  ifstream f_e(filename4);
  ifstream f_nb(filename5);
        
  if (!f_e.good()) {
    cout << "I/O error with " << filename4 << endl;
    exit(1);
  }
  if (!f_nb.good()){
    cout << "I/O error with " << filename5 << endl;
    exit(1);
  }
        
  nT     = 771;
  nmuB   = 451;
  tmin   = 0.030;
  tmax   = 0.800;
  mubmax = 0.450;
  mubmin = 0.0;
        
  epsilontab = new double* [nT];
  nbtab      = new double* [nT];
  for (int i = 0; i < nT; i++) {
    epsilontab[i] = new double[nmuB];
    nbtab[i]      = new double[nmuB];
  }
  double temp;
  for(int it = 0; it < nT; it++)
     for(int imub = 0; imub < nmuB; imub++){
        f_e   >> temp >> temp >> epsilontab[it][imub];
        f_nb  >> temp >> temp >> nbtab[it][imub];
        epsilontab[it][imub] = epsilontab[it][imub] * pow((it+30)*0.001, 4) / pow(hbarC, 3); //-> GeV/fm^3
        nbtab[it][imub]      = nbtab[it][imub]      * pow((it+30)*0.001, 3) / pow(hbarC, 3); //-> 1/fm^3
    }
}
EoSauxBEST::EoSauxBEST(const EoSauxBEST& ee){
    emax    = ee.emax;
    nmax    = ee.nmax;
    emin    = ee.emin;
    nmin    = ee.nmin;
    ne      = ee.ne;
    nn      = ee.nn;
    ptab    = new double* [ne];
    Ttab    = new double* [ne];
    mubtab  = new double* [ne];
    mustab  = new double* [ne];
    stab    = new double* [ne];
    
    for (int i = 0; i < ne; i++) {
        ptab[i] = new double[nn];
        Ttab[i] = new double[nn];
        mubtab[i] = new double[nn];
        mustab[i] = new double[nn];
        stab[i] = new double[nn];
    }
    for (int i = 0; i < ne; i++) {
        for (int j = 0; j < nn; j++) {
            ptab[i][j] = ee.ptab[i][j];
            Ttab[i][j] = ee.Ttab[i][j];
            mubtab[i][j] = ee.mubtab[i][j];
            stab[i][j] = ee.stab[i][j];
            mustab[i][j] = ee.mustab[i][j];
    }}
 
    nT = ee.nT;
    nmuB = ee.nmuB;
    tmin   = ee.tmin;
    tmax   = ee.tmax;
    mubmax = ee.mubmax;
    mubmin = ee.mubmin;
    epsilontab = new double* [nT];
    nbtab = new double* [nT];
    for (int i = 0; i < nT; i++) {
      epsilontab[i] = new double[nmuB];
      nbtab[i] = new double[nmuB];
    }
    for( int it = 0; it < nT; it++)
        for(int imub = 0; imub < nmuB; imub++){
            epsilontab[it][imub] = ee.epsilontab[it][imub];
            nbtab[it][imub] = ee.nbtab[it][imub];
    }
}
EoSauxBEST::~EoSauxBEST() {
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

  for (int i = 0; i < nT; i++) {
    delete[] epsilontab[i];
    delete[] nbtab[i];
  }
  delete[] epsilontab;
  delete[] nbtab;
}

void EoSauxBEST::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {
  if (e < 0.) {
    T = mub = mus = p = 0.0;
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
  const double em = e - emin - ie * de;
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
}

double EoSauxBEST::p(double e, double nb) {
  if (e < 1e-10) return 0.0;
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
  for (int je = 0; je < 2; je++)
    for (int jn = 0; jn < 2; jn++)
        p += we[je] * wn[jn] * ptab[ie + je][in + jn];
 // if (e < 1e-15 && nb < 1e-15 ) { p = ptab[0][0];}
//  if (p < 1e-15) p = 0.0;
  return p;
}

void EoSauxBEST::geteosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) {
  
    e = 0.0; n_b= 0.0; n_q= 0.0; n_s= 0.0; p = 0.0;
    double epsilon = 0.0; double nb = 0.0;
    if(T < 0) return;
    const double dt = (tmax - tmin) / (nT - 1);
    const double dmub = (mubmax - mubmin) / (nmuB - 1);
    int it = (int)((T - tmin) / dt);
    int imub = (int)((mu_b - mubmin) / dmub);
    if (it < 0) it = 0;
    if (imub < 0) imub = 0;
    if (it > nT - 2) it = nT - 2;
    if (imub > nmuB - 2) imub = nmuB - 2;
    const double tm   = T - tmin - it * dt;
    const double mubm = mu_b - mubmin - imub * dmub;
    double wt[2] = {1. - tm / dt, tm / dt};
    double wm[2] = {1. - mubm / dmub, mubm / dmub};

    
    if(imub == 0){ // on the "left" edge of the table
        for (int je = 0; je < 2; je++){
            epsilon += wt[je]  * epsilontab[it + je][0];
            nb += wt[je]  * nbtab[it + je][0];
        }
    }
    else if(it == 0){
        for (int jn = 0; jn < 2; jn++)
        {
            epsilon += wm[jn]  * epsilontab[0][imub + jn];
            nb      += wm[jn]  * nbtab[0][imub + jn];
        }
    }
    else{
        for (int je = 0; je < 2; je++)
            for (int jn = 0; jn < 2; jn++){
                      epsilon += wt[je] * wm[jn] * epsilontab[it + je][imub + jn];
                      nb      += wt[je] * wm[jn] * nbtab[it + je][imub + jn];
    }}
    e   = epsilon;
    n_b = nb;
}
/*
double EoSauxBEST::get_nb(double T, double mub){

    if(T < 0) return 0.0;
    
    const double dt = (tmax - tmin) / (nT - 1);
    const double dmub = (mubmax - mubmin) / (nmuB - 1);
    int it = (int)((T - tmin) / dt);
    int imub = (int)((mub - mubmin) / dmub);
    if (it < 0)   it   = 0;
    if (imub < 0) imub = 0;
    if (it > nT - 2) it = nT - 2;
    if (imub > nmuB - 2) imub = nmuB - 2;
    const double tm = T - tmin - it * dt;
    const double mubm = mub - mubmin - imub * dmub;
    double wt[2] = {1. - tm / dt, tm / dt};
    double wm[2] = {1. - mubm / dmub, mubm / dmub};

    double nb = 0.0;
    if(imub == 0){ // on the "left" edge of the table
      for (int je = 0; je < 2; je++)
             
    }else if(it == 0){
        for (int jn = 0; jn < 2; jn++)
            nb += wm[jn]  * nbtab[0][imub + jn];
    }else{
      for (int je = 0; je < 2; je++)
          for (int jn = 0; jn < 2; jn++)
               nb += wt[je] * wm[jn] * nbtab[it + je][imub + jn];
    }
    return nb;
}

*/
eoBEST::eoBEST(char* filename) {
    for(int i = 0; i < 6; i++){
        EoSauxBEST object(filename, i);
        eos_t.push_back(object);
    }
}

eoBEST::~eoBEST() {
}

void eoBEST::eosranges(double &_emax, double &_e0, double &_nmax,
                          double &_n0, int &_ne, int &_nn) {
 _emax = eos_t[5].get_emax();
 _nmax = eos_t[5].get_nmax();
 _e0 = 0.;  // I don't think it is a meaningul param, so return 0
 _n0 = 0.;  // I don't think it is a meaningul param, so return 0
 _ne = eos_t[5].getne(); // This does not make sense either since there are two EoS
 _nn = eos_t[5].getnn(); // tables, so I return the size of larger table.
}

void eoBEST::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
    int iitab = -99;
    if(e <= 1e-12){
        T = 0.0;
        p = 0.0;
        e = 0.0;
    }
    for(int i = 0; i <6; i++){
        if (e >= ed_bounds[i] && e < ed_bounds[i+1] && nb >= 0.0 && nb < nB_bounds[i+1])  {
            eos_t[i].get(e, nb, p, T, mub, mus);
            iitab = i;
            break;
        }
    }
    if (iitab == -99){
        p = 0.2964 * e;
        T = 0.15120476935 * pow(e, 0.25);
        mub = mus = 0.0;
    }
  muq = 0.0;
}

double eoBEST::p(double e, double nb, double nq, double ns) {
  if(e <= 1e-12) return 0.0;
  for(int i = 0; i < 6; i++){
      if (e >= ed_bounds[i] && e < ed_bounds[i+1] && nb >= 0.0 && nb < nB_bounds[i+1])  {
          return eos_t[i].p(e, nb);
      }
  }
  return 0.2964 * e;
}

void eoBEST::eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) {
    eos_t[0].geteosorginal(T, mu_b, mu_q, mu_s, e, n_b, n_q, n_s, p);
    
}
