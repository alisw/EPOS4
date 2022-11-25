//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

//#include <TError.h>
//#include <TApplication.h>
//#include <TGraph.h>
//#include <TCanvas.h>
//#include <TMath.h>
//#include <TGraph.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
//#include <TF1.h>

#include "eos.h"
#include "eoCEM.h"

using namespace std;

// ---- auxiliary EoS class. Two instances (objects) of this class will be
// created
// ---- to store two EoS tables: chiraleos and chiralsmall
// ---- Not performance-wise, but the code is elegant :)
class EoSauxCEM {
  double emax, nmax, emin, nmin;
  int ne, nn; //number of epsilon and nb points
  double** ptab, **Ttab, **mubtab, **mustab, **stab;

 public:
  EoSauxCEM(char* filename, int Ne, int Nn);
  ~EoSauxCEM();
  double get_emax() { return emax; }
  double get_nmax() { return nmax; }
  double getne() { return ne; }
  double getnn() { return nn; }
  void get(double e, double nb, double& p, double& T, double& mub, double& mus);
  double p(double e, double nb);
};

EoSauxCEM::EoSauxCEM(char* filename, int Ne, int Nn) {
  ne = Ne;
  nn = Nn;
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
  ifstream fin(filename);
  double a;  // dummy variable
  if (!fin.good()) {
    cout << "I/O error with " << filename << endl;
    exit(1);
  }
  for (int in = 0; in < nn; in++)
    for (int ie = 0; ie < ne; ie++) {
      // T [MeV], mu_q [MeV], energy density [e_0], pressure [e_0], baryon
      // density [n_0], entropy density [n_0], mu_S [MeV], placeholder (not
      // important).
      fin >> Ttab[ie][in] >> mubtab[ie][in] >> e[ie] >> ptab[ie][in] >> n[in] >>
          stab[ie][in] >> mustab[ie][in] >> a;
      Ttab[ie][in] /= 1000.0;    // --> T[GeV]
      mubtab[ie][in] /= 1000.0;  // --> mub[GeV]
      mustab[ie][in] /= 1000.0;  // --> mus[GeV]
      ptab[ie][in] *= 0.146;     // --> p[GeV/fm3]
      stab[ie][in] *= 0.15;      // --> s[1/fm3]
    }
    ///////---------------------
  emin = e[0] * 0.146;
  emax = e[ne - 1] * 0.146;
  nmin = n[0] * 0.15;
  nmax = n[nn - 1] * 0.15; // 
  /////////////////
  cout << "EoSaux: table " << filename
       << " read, [emin,emax,nmin,nmax] = " << emin << "  " << emax << "  "
       << nmin << "  " << nmax << endl;
  delete[] e;
  delete[] n;
}

EoSauxCEM::~EoSauxCEM() {
  for (int i = 0; i < ne; i++) {
    delete[] ptab[i];
    delete[] Ttab[i];
    delete[] mubtab[i];
    delete[] mustab[i];
    delete[] stab[i];
  }
  delete ptab;
  delete Ttab;
  delete mubtab;
  delete mustab;
  delete stab;
}

void EoSauxCEM::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {
  if (e < 0.) {
    T = mub = mus = p = 0.;
    return;
  }
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

  T = mub = mus = p = 0.0;
  for (int je = 0; je < 2; je++)
    for (int jn = 0; jn < 2; jn++) {
      p += we[je] * wn[jn] * ptab[ie + je][in + jn];
      T += we[je] * wn[jn] * Ttab[ie + je][in + jn];
      mub += we[je] * wn[jn] * mubtab[ie + je][in + jn];
      mus += we[je] * wn[jn] * mustab[ie + je][in + jn];
    }
  if (p < 0.0) p = 0.0;
  // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}

double EoSauxCEM::p(double e, double nb) {
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
  for (int je = 0; je < 2; je++)
    for (int jn = 0; jn < 2; jn++)
      p += we[je] * wn[jn] * ptab[ie + je][in + jn];

  if (p < 0.0) p = 0.0;
  return p;
  // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}

eoCEM::eoCEM(char* fileMini, char* fileSmall, char* fileLarge) {
  eosmini  = new EoSauxCEM(fileMini, 201, 201);
  eossmall = new EoSauxCEM(fileSmall, 201, 201);
  eosbig   = new EoSauxCEM(fileLarge, 2001, 401);
}

eoCEM::~eoCEM() {
  delete eosbig;
  delete eossmall;
  delete eosmini;
}

void eoCEM::eosranges(double &_emax, double &_e0, double &_nmax,
                          double &_n0, int &_ne, int &_nn) {
 _emax = eosbig->get_emax();
 _nmax = eosbig->get_nmax();
 _e0 = 0.;  // I don't think it is a meaningul param, so return 0
 _n0 = 0.;  // I don't think it is a meaningul param, so return 0
 _ne = eosbig->getne(); // This does not make sense either since there are two EoS
 _nn = eosbig->getnn(); // tables, so I return the size of larger table.
}

void eoCEM::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
 /* if (e < 0.1*0.146 && nb < 0.02*0.159)
    eosmini->get(e, nb, p, T, mub, mus);
  else if (e < 10*0.146 && nb < 2.0*0.159)
    eossmall->get(e, nb, p, T, mub, mus);
  else if (e < 1000*0.146 && nb < 40.0*0.159)
    eosbig->get(e, nb, p, T, mub, mus);
  else {*/
   if (e < 0.1 && nb < 0.02)
    eosmini->get(e, nb, p, T, mub, mus);
  else if (e < 10 && nb < 2.0)
    eossmall->get(e, nb, p, T, mub, mus);
  else if (e < 1000&& nb < 40.0)
    eosbig->get(e, nb, p, T, mub, mus);
  else {
    p = 0.2964 * e;
    T = 0.15120476935 * pow(e, 0.25);
    mub = mus = 0.0;
  }
  muq = 0.0;  // generally it's not zero, but...but
}

double eoCEM::p(double e, double nb, double nq, double ns) {
  if (e < 146. && nb < 6.)
    return eosbig->p(e, nb);
  else
    return 0.2964 * e;
}

void eoCEM::eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) {
    cout << " >>>>  no original table avaliable <<<< " << endl;
}
