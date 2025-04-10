#include "Eosp.hpp"

//! default class constructor 
Eosp::Eosp() {}

//! default class destructor 
Eosp::~Eosp() {
}

// map the fortran common block ceost2
extern "C" {
extern struct {
  float uEeos;
  float oEeos;
} ceost2_;
}

// get value from fortran common block
void Eosp::getFortranOEeos() { this->oEeos = ceost2_.oEeos; }

// set value from fortran common block
void Eosp::setFortranOEeos() { ceost2_.oEeos = this->oEeos; }

// c++ attribute accessor methods
//! get oEeos attribute value
/*!
\return The oEeos attribute value - float type
*/
float Eosp::getOEeos() {
  // get value from fortran common block
  this->getFortranOEeos();
  return this->oEeos;
}

// c++ attribute mutator methods
//! set oEeos attribute value
/*!
\param oEeosParam - the oEeos attribute value - float type
*/
void Eosp::setOEeos(float oEeosParam) {
  this->oEeos = oEeosParam;
  // set value from fortran common block
  this->setFortranOEeos();
}

// Eosp instance pointer
Eosp *eosp;

// Fortran / C wrapper functions
extern "C" {
//! create Eosp instance
void createeosp_() { eosp = new Eosp(); }

//! delete Eosp instance
void destroyeosp_() { delete eosp; }

//! get oEeos attribute value
/*!
\return The oEeos attribute value - float type
*/
float geteospoeeos_() { return eosp->getOEeos(); }

//! set oEeos attribute value
/*!
\param oEeosParam - the oEeos attribute value - float type
*/
void seteospoeeos_(float *oEeosParam) { eosp->setOEeos(*oEeosParam); }

}
