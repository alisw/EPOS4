//
// This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
// (See COPYING file for the text of the licence)
//

#include <iostream>

extern "C" {
void getnptl_  (int * );
void getiorptl_(int *, int *);
void getjorptl_(int *, int *);
void getifrptl_(int *, int *, int *);
void getidptl_ (int *, int *);
void getistptl_(int *, int *);
void getityptl_(int *, int *);
void getcharge_(int *, float *);
void getihadron_(int *, int *);
void getpptl_  (int *, float *, float *, float *, float *, float *);
void getxorptl_(int *, float *, float *, float *, float *);
void getbim_(float *bim);
void idepos2pdg_(int *, int *);
float polar_(float *x,float *y);
}


/**
 * Get number of particles for current event
 *
 */
void getNptl(int *nptl) { 
  getnptl_(nptl); 
}
/**
 * Get ior for particle i
 *
 */
void getIorptl(int *i, int *ior) {
  getiorptl_(i, ior);
}
/**
 * Get jor for particle i
 *
 */
void getJorptl(int *i, int *jor) {
  getjorptl_(i, jor);
}
/**
 * Get ifr1, ifr2 for particle i
 *
 */
void getIfrptl(int *i, int *ifr1, int *ifr2) {
  getifrptl_(i, ifr1, ifr2);
}
/**
 * Get id for particle i
 *
 */
void getIdptl(int *i, int *id) {
  getidptl_(i, id);
}
/**
 * Get ist for particle i
 *
 */
void getIstptl(int *i, int *ist) {
  getistptl_(i, ist);
}
/**
 * Get ity for particle i
 *
 */
void getItyptl(int *i, int *ity) {
  getityptl_(i, ity);
}

/**
 * Get charge for particle i
 *
 */
void getCharge(int *i, float *charge) {
  getcharge_(i, charge);
}
/**
 * Get ihadron for particle i
 *
 */
void getIhadron(int *i, int *ihadron) {
  getihadron_(i, ihadron);
}
/**
 * Get 3-momentum, energy, mass for particle i
 * @param i particle index
 * @param px,py,pz -- 3-momentum of particle
 * @param energy -- energy of particle
 * @param mass -- mass of particle
 */
void getPptl(int *i, float *px, float *py, float *pz, float *energy, float *mass) {
  getpptl_(i,px,py,pz,energy,mass);
}
/**
 * Get formation position and time for particle i
 * @param i particle index
 * @param x,y,z -- position of particle
 * @param t -- formation time of particle
 */
void getXorptl(int *i, float *x, float *y, float *z, float *t) {
  getxorptl_(i,x,y,z,t);
}

/**
 * Get impact parameter for current event
 * @param bim impact parameter
 */
void getBim(float *bim) { 
  getbim_(bim);
}

/**
 * Get PDG id code from the EPOS one
 * @param idepos EPOS id code
 * @param idpdg  PDG id code
 */
void getIdpdg(int *idepos, int *idpdg) { 
  idepos2pdg_(idepos,idpdg);
}
/**
 * Get polar angle in rad for given x and y
 *
 */
void getPolar(float *x, float *y, float *phi) {
  *phi=polar_(x,y);
}

