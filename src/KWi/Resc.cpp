#include "Resc.hpp"

//! default class constructor 
Resc::Resc() {}

//! default class destructor 
Resc::~Resc() {
}

// map the fortran common block resc1
extern "C" {
extern struct {
  float tauZero;
  float tauOne;
  float tauTwo;
  float deltau;
  float numtau;
  float amsiac;
  float amprif;
  float tauUp;
  float tauThree;
} resc1_;
}

// get value from fortran common block
void Resc::getFortranTauZero() { this->tauZero = resc1_.tauZero; }

void Resc::getFortranTauOne() { this->tauOne = resc1_.tauOne; }

void Resc::getFortranTauTwo() { this->tauTwo = resc1_.tauTwo; }

void Resc::getFortranTauUp() { this->tauUp = resc1_.tauUp; }

void Resc::getFortranTauThree() { this->tauThree = resc1_.tauThree; }

// set value from fortran common block
void Resc::setFortranTauZero() { resc1_.tauZero = this->tauZero; }

void Resc::setFortranTauOne() { resc1_.tauOne = this->tauOne; }

void Resc::setFortranTauTwo() { resc1_.tauTwo = this->tauTwo; }

void Resc::setFortranTauUp() { resc1_.tauUp = this->tauUp; }

void Resc::setFortranTauThree() { resc1_.tauThree = this->tauThree; }

// map the fortran common block resc3
extern "C" {
extern struct {
  float dscale;
  float cepara;
  int iceopt;
  float delamf;
  float deuamf;
  float etaos;
  float zetaos;
} resc3_;
}

// get value from fortran common block
void Resc::getFortranEtaos() { this->etaos = resc3_.etaos; }

void Resc::getFortranZetaos() { this->zetaos = resc3_.zetaos; }

// set value from fortran common block
void Resc::setFortranEtaos() { resc3_.etaos = this->etaos; }

void Resc::setFortranZetaos() { resc3_.zetaos = this->zetaos; }

// c++ attribute accessor methods
//! get tauZero attribute value
/*!
\return The tauZero attribute value - float type
*/
float Resc::getTauZero() {
  // get value from fortran common block
  this->getFortranTauZero();
  return this->tauZero;
}

//! get tauOne attribute value
/*!
\return The tauOne attribute value - float type
*/
float Resc::getTauOne() {
  // get value from fortran common block
  this->getFortranTauOne();
  return this->tauOne;
}

//! get tauTwo attribute value
/*!
\return The tauTwo attribute value - float type
*/
float Resc::getTauTwo() {
  // get value from fortran common block
  this->getFortranTauTwo();
  return this->tauTwo;
}

//! get tauThree attribute value
/*!
\return The tauThree attribute value - float type
*/
float Resc::getTauThree() {
  // get value from fortran common block
  this->getFortranTauThree();
  return this->tauThree;
}

//! get tauUp attribute value
/*!
\return The tauUp attribute value - float type
*/
float Resc::getTauUp() {
  // get value from fortran common block
  this->getFortranTauUp();
  return this->tauUp;
}

//! get etaos attribute value
/*!
\return The etaos attribute value - float type
*/
float Resc::getEtaos() {
  // get value from fortran common block
  this->getFortranEtaos();
  return this->etaos;
}

//! get zetaos attribute value
/*!
\return The zetaos attribute value - float type
*/
float Resc::getZetaos() {
  // get value from fortran common block
  this->getFortranZetaos();
  return this->zetaos;
}

// c++ attribute mutator methods
//! set tauZero attribute value
/*!
\param tauZeroParam - the tauZero attribute value - float type
*/
void Resc::setTauZero(float tauZeroParam) {
  this->tauZero = tauZeroParam;
  // set value from fortran common block
  this->setFortranTauZero();
}

//! set tauOne attribute value
/*!
\param tauOneParam - the tauOne attribute value - float type
*/
void Resc::setTauOne(float tauOneParam) {
  this->tauOne = tauOneParam;
  // set value from fortran common block
  this->setFortranTauOne();
}

//! set tauTwo attribute value
/*!
\param tauTwoParam - the tauTwo attribute value - float type
*/
void Resc::setTauTwo(float tauTwoParam) {
  this->tauTwo = tauTwoParam;
  // set value from fortran common block
  this->setFortranTauTwo();
}

//! set tauThree attribute value
/*!
\param tauThreeParam - the tauThree attribute value - float type
*/
void Resc::setTauThree(float tauThreeParam) {
  this->tauThree = tauThreeParam;
  // set value from fortran common block
  this->setFortranTauThree();
}

//! set tauUp attribute value
/*!
\param tauUpParam - the tauUp attribute value - float type
*/
void Resc::setTauUp(float tauUpParam) {
  this->tauUp = tauUpParam;
  // set value from fortran common block
  this->setFortranTauUp();
}

//! set etaos attribute value
/*!
\param etaosParam - the etaos attribute value - float type
*/
void Resc::setEtaos(float etaosParam) {
  this->etaos = etaosParam;
  // set value from fortran common block
  this->setFortranEtaos();
}

//! set zetaos attribute value
/*!
\param zetaosParam - the zetaos attribute value - float type
*/
void Resc::setZetaos(float zetaosParam) {
  this->zetaos = zetaosParam;
  // set value from fortran common block
  this->setFortranZetaos();
}

//! get tauZero, tauOne, tauTwo, tauThree, tauUp attribute values
/*!
\param tauZeroParam - the tauZero attribute address - float pointer type
\param tauOneParam - the tauOne attribute address - float pointer type
\param tauTwoParam - the tauTwo attribute address - float pointer type
\param tauThreeParam - the tauThree attribute address - float pointer type
\param tauUpParam - the tauUp attribute address - float pointer type
*/
void Resc::getTau(float *tauZeroParam, float *tauOneParam, float *tauTwoParam, float *tauThreeParam, float *tauUpParam) {
  *tauZeroParam = this->getTauZero();
  *tauOneParam = this->getTauOne();
  *tauTwoParam = this->getTauTwo();
  *tauThreeParam = this->getTauThree();
  *tauUpParam = this->getTauUp();
}

//! set tauZero, tauOne, tauTwo, tauThree, tauUp attribute values
/*!
\param tauZeroParam - the tauZero attribute address - float type
\param tauOneParam - the tauOne attribute address - float type
\param tauTwoParam - the tauTwo attribute address - float type
\param tauThreeParam - the tauThree attribute address - float type
\param tauUpParam - the tauUp attribute address - float type
*/
void Resc::setTau(float tauZeroParam, float tauOneParam, float tauTwoParam, float tauThreeParam, float tauUpParam) {
  this->setTauZero(tauZeroParam);
  this->setTauOne(tauOneParam);
  this->setTauTwo(tauTwoParam);
  this->setTauThree(tauThreeParam);
  this->setTauUp(tauUpParam);
}

// Resc instance pointer
Resc *resc;

// Fortran / C wrapper functions
extern "C" {
//! create Resc instance
void createresc_() { resc = new Resc(); }

//! delete Resc instance
void destroyresc_() { delete resc; }

//! get tauZero attribute value
/*!
\return The tauZero attribute value - float type
*/
float getresctauzero_() { return resc->getTauZero(); }

//! set tauZero attribute value
/*!
\param tauZeroParam - the tauZero attribute value - float type
*/
void setresctauzero_(float *tauZeroParam) { resc->setTauZero(*tauZeroParam); }

//! get tauOne attribute value
/*!
\return The tauOne attribute value - float type
*/
float getresctauone_() { return resc->getTauOne(); }

//! set tauOne attribute value
/*!
\param tauOneParam - the tauOne attribute value - float type
*/
void setresctauone_(float *tauOneParam) { resc->setTauOne(*tauOneParam); }

//! get tauTwo attribute value
/*!
\return The tauTwo attribute value - float type
*/
float getresctautwo_() { return resc->getTauTwo(); }

//! set tauTwo attribute value
/*!
\param tauTwoParam - the tauTwo attribute value - float type
*/
void setresctautwo_(float *tauTwoParam) { resc->setTauTwo(*tauTwoParam); }

//! get tauThree attribute value
/*!
\return The tauThree attribute value - float type
*/
float getresctauthree_() { return resc->getTauThree(); }

//! set tauThree attribute value
/*!
\param tauThreeParam - the tauThree attribute value - float type
*/
void setresctauthree_(float *tauThreeParam) { resc->setTauThree(*tauThreeParam); }

//! get tauUp attribute value
/*!
\return The tauUp attribute value - float type
*/
float getresctauup_() { return resc->getTauUp(); }

//! set tauUp attribute value
/*!
\param tauUpParam - the tauUp attribute value - float type
*/
void setresctauup_(float *tauUpParam) { resc->setTauUp(*tauUpParam); }

//! get etaos attribute value
/*!
\return The etaos attribute value - float type
*/
float getrescetaos_() { return resc->getEtaos(); }

//! set etaos attribute value
/*!
\param etaosParam - the etaos attribute value - float type
*/
void setrescetaos_(float *etaosParam) { resc->setEtaos(*etaosParam); }

//! get zetaos attribute value
/*!
\return The zetaos attribute value - float type
*/
float getresczetaos_() { return resc->getZetaos(); }

//! set zetaos attribute value
/*!
\param zetaosParam - the zetaos attribute value - float type
*/
void setresczetaos_(float *zetaosParam) { resc->setZetaos(*zetaosParam); }

//! get tauZero, tauOne, tauTwo, tauThree, tauUp attribute values
/*!
\param tauZeroParam - the tauZero attribute address - float pointer type
\param tauOneParam - the tauOne attribute address - float pointer type
\param tauTwoParam - the tauTwo attribute address - float pointer type
\param tauThreeParam - the tauThree attribute address - float pointer type
\param tauUpParam - the tauUp attribute address - float pointer type
*/
void getresctau_(float* tauZeroParam, float* tauOneParam, float* tauTwoParam, float* tauThreeParam, float* tauUpParam) {
  resc->getTau(tauZeroParam,tauOneParam,tauTwoParam,tauThreeParam,tauUpParam);
}

//! set tauZero, tauOne, tauTwo, tauThree, tauUp attribute values
/*!
\param tauZeroParam - the tauZero attribute address - float pointer type
\param tauOneParam - the tauOne attribute address - float pointer type
\param tauTwoParam - the tauTwo attribute address - float pointer type
\param tauThreeParam - the tauThree attribute address - float pointer type
\param tauUpParam - the tauUp attribute address - float pointer type
*/
void setresctau_(float* tauZeroParam, float* tauOneParam, float* tauTwoParam, float* tauThreeParam, float* tauUpParam) {
  resc->setTau(*tauZeroParam,*tauOneParam,*tauTwoParam,*tauThreeParam,*tauUpParam);
}

}
