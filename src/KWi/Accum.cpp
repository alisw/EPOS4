#include "Accum.hpp"

//! default class constructor 
Accum::Accum() {}

//! default class destructor 
Accum::~Accum() {
}

// map the fortran common block accum
extern "C" {
extern struct {
  int imsg;
  int ntevt;
  int nrevt;
} accum_;
}

// get value from fortran common block
void Accum::getFortranImsg() { this->imsg = accum_.imsg; }

void Accum::getFortranNtevt() { this->ntevt = accum_.ntevt; }

void Accum::getFortranNrevt() { this->nrevt = accum_.nrevt; }

// set value from fortran common block
void Accum::setFortranImsg() { accum_.imsg = this->imsg; }

void Accum::setFortranNtevt() { accum_.ntevt = this->ntevt; }

void Accum::setFortranNrevt() { accum_.nrevt = this->nrevt; }

// c++ attribute accessor methods
//! get imsg attribute value
/*!
\return The imsg attribute value - int type
*/
int Accum::getImsg() {
  // get value from fortran common block
  this->getFortranImsg();
  return this->imsg;
}

//! get jerr attribute value
/*!
\return The jerr attribute value - Array<int> * type
*/
Array<int> *Accum::getJerr() { return this->jerr; }

//! get ntevt attribute value
/*!
\return The ntevt attribute value - int type
*/
int Accum::getNtevt() {
  // get value from fortran common block
  this->getFortranNtevt();
  return this->ntevt;
}

//! get nrevt attribute value
/*!
\return The nrevt attribute value - int type
*/
int Accum::getNrevt() {
  // get value from fortran common block
  this->getFortranNrevt();
  return this->nrevt;
}

// c++ attribute mutator methods
//! set imsg attribute value
/*!
\param imsgParam - the imsg attribute value - int type
*/
void Accum::setImsg(int imsgParam) {
  this->imsg = imsgParam;
  // set value from fortran common block
  this->setFortranImsg();
}

//! set jerr attribute value
/*!
\param jerrParam - the jerr attribute value - Array<int> * type
*/
void Accum::setJerr(Array<int> *jerrParam) { this->jerr = jerrParam; }

//! set ntevt attribute value
/*!
\param ntevtParam - the ntevt attribute value - int type
*/
void Accum::setNtevt(int ntevtParam) {
  this->ntevt = ntevtParam;
  // set value from fortran common block
  this->setFortranNtevt();
}

//! set nrevt attribute value
/*!
\param nrevtParam - the nrevt attribute value - int type
*/
void Accum::setNrevt(int nrevtParam) {
  this->nrevt = nrevtParam;
  // set value from fortran common block
  this->setFortranNrevt();
}

// c++ array attribute methods
//! create jerr attribute array of int type
/*!
\param i - indice - int type
*/
void Accum::createJerr(int i) {
  int dim = 1;
  int shape[1] = {i};
  this->jerr = new Array<int>(dim, shape);
}

//! get jerr attribute array value
/*!
\param i - indice - int type
\return The jerr attribute value - int type
*/
int Accum::getJerrValue(int i){
  int index[1] = {i};
  return this->jerr->get(index);
}

//! set jerr attribute array value
/*!
\param i - indice - int type
\param value - the value to be set - int type
*/
void Accum::setJerrValue(int i, int value) {
  int index[1] = {i};
  this->jerr->set(index, value);
}

//! create jerr initialize method
void Accum::initializeJerr(int value) {
  this->jerr->initialize(value);
}

// Accum instance pointer
Accum *accum;

// Fortran / C wrapper functions
extern "C" {
//! create Accum instance
void createaccum_() { accum = new Accum(); }

//! delete Accum instance
void destroyaccum_() { delete accum; }

//! get imsg attribute value
/*!
\return The imsg attribute value - int type
*/
int getaccumimsg_() { return accum->getImsg(); }

//! set imsg attribute value
/*!
\param imsgParam - the imsg attribute value - int type
*/
void setaccumimsg_(int *imsgParam) { accum->setImsg(*imsgParam); }

//! create jerr array attribute
/*!\param i - indice - int type
*/void createaccumjerr_(int* i) {
  accum->createJerr(*i);
}

//! initialize jerr array attribute
void initializeaccumjerr_(int *value) {
  accum->initializeJerr(*value);
}

//! set jerr array value
/*!
\param i - indice - int pointer type
\param value - value to set - int pointer type
*/
void setaccumjerr_(int *i, int *value) {
  accum->setJerrValue(*i, *value);
}

//! get jerr array value
/*!
\param i indice - int pointer type
\return value - value to get - int pointer type
*/
int getaccumjerr_(int *i) {
  return accum->getJerrValue(*i);
}

//! delete jerr array
void destroyaccumjerr_() { delete accum->getJerr(); }

//! get ntevt attribute value
/*!
\return The ntevt attribute value - int type
*/
int getaccumntevt_() { return accum->getNtevt(); }

//! set ntevt attribute value
/*!
\param ntevtParam - the ntevt attribute value - int type
*/
void setaccumntevt_(int *ntevtParam) { accum->setNtevt(*ntevtParam); }

//! get nrevt attribute value
/*!
\return The nrevt attribute value - int type
*/
int getaccumnrevt_() { return accum->getNrevt(); }

//! set nrevt attribute value
/*!
\param nrevtParam - the nrevt attribute value - int type
*/
void setaccumnrevt_(int *nrevtParam) { accum->setNrevt(*nrevtParam); }

}
