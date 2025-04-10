#include "Core.hpp"

//! default class constructor 
Core::Core() {}

//! default class destructor 
Core::~Core() {
}

// map the fortran common block core1
extern "C" {
extern struct {
  int icotabm;
  int icotabr;
  int icocore;
  int jcorona;
} core1_;
}

// get value from fortran common block
void Core::getFortranIcotabr() { this->icotabr = core1_.icotabr; }

// set value from fortran common block
void Core::setFortranIcotabr() { core1_.icotabr = this->icotabr; }

// c++ attribute accessor methods
//! get icotabr attribute value
/*!
\return The icotabr attribute value - int type
*/
int Core::getIcotabr() {
  // get value from fortran common block
  this->getFortranIcotabr();
  return this->icotabr;
}

// c++ attribute mutator methods
//! set icotabr attribute value
/*!
\param icotabrParam - the icotabr attribute value - int type
*/
void Core::setIcotabr(int icotabrParam) {
  this->icotabr = icotabrParam;
  // set value from fortran common block
  this->setFortranIcotabr();
}

// Core instance pointer
Core *core;

// Fortran / C wrapper functions
extern "C" {
//! create Core instance
void createcore_() { core = new Core(); }

//! delete Core instance
void destroycore_() { delete core; }

//! get icotabr attribute value
/*!
\return The icotabr attribute value - int type
*/
int getcoreicotabr_() { return core->getIcotabr(); }

//! set icotabr attribute value
/*!
\param icotabrParam - the icotabr attribute value - int type
*/
void setcoreicotabr_(int *icotabrParam) { core->setIcotabr(*icotabrParam); }

}
