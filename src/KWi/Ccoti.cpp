#include "Ccoti.hpp"

//! default class constructor 
Ccoti::Ccoti() {}

//! default class destructor 
Ccoti::~Ccoti() {
}

// c++ attribute accessor methods
//! get coti attribute value
/*!
\return The coti attribute value - Array<float> * type
*/
Array<float> *Ccoti::getCoti() { return this->coti; }

// c++ attribute mutator methods
//! set coti attribute value
/*!
\param cotiParam - the coti attribute value - Array<float> * type
*/
void Ccoti::setCoti(Array<float> *cotiParam) { this->coti = cotiParam; }

// c++ array attribute methods
//! create coti attribute array of float type
/*!
\param x - indice - int type
\param y - indice - int type
*/
void Ccoti::createCoti(int x, int y) {
  int dim = 2;
  int shape[2] = {x, y};
  this->coti = new Array<float>(dim, shape);
}

//! get coti attribute array value
/*!
\param x - indice - int type
\param y - indice - int type
\return The coti attribute value - float type
*/
float Ccoti::getCotiValue(int x, int y){
  int index[2] = {x, y};
  return this->coti->get(index);
}

//! set coti attribute array value
/*!
\param x - indice - int type
\param y - indice - int type
\param value - the value to be set - float type
*/
void Ccoti::setCotiValue(int x, int y, float value) {
  int index[2] = {x, y};
  this->coti->set(index, value);
}

//! create coti initialize method
void Ccoti::initializeCoti(float value) {
  this->coti->initialize(value);
}

// Ccoti instance pointer
Ccoti *ccoti;

// Fortran / C wrapper functions
extern "C" {
//! create Ccoti instance
void createccoti_() { ccoti = new Ccoti(); }

//! delete Ccoti instance
void destroyccoti_() { delete ccoti; }

//! create coti array attribute
/*!\param x - indice - int type
\param y - indice - int type
*/void createccoticoti_(int* x, int* y) {
  ccoti->createCoti(*x, *y);
}

//! initialize coti array attribute
void initializeccoticoti_(float *value) {
  ccoti->initializeCoti(*value);
}

//! set coti array value
/*!
\param x - indice - int pointer type
\param y - indice - int pointer type
\param value - value to set - float pointer type
*/
void setccoticoti_(int *x, int *y, float *value) {
  ccoti->setCotiValue(*x, *y, *value);
}

//! get coti array value
/*!
\param x indice - int pointer type
\param y indice - int pointer type
\return value - value to get - float pointer type
*/
float getccoticoti_(int *x, int *y) {
  return ccoti->getCotiValue(*x, *y);
}

//! delete coti array
void destroyccoticoti_() { delete ccoti->getCoti(); }

}
