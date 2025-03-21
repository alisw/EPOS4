#include "Outlist.hpp"

//! default class constructor 
Outlist::Outlist() {}

//! default class destructor 
Outlist::~Outlist() {
}

// c++ array attribute methods
//! create thename attribute array of char type
/*!
\param i - indice - int type
*/
void Outlist::createThename(int i) {
  int dim = 1;
  int shape[1] = {i};
  this->thename = new Array<char>(dim, shape);
}

//########################################### get/set #############################################

// c++ attribute accessor methods
//! get toscreen attribute value
/*!
\return The toscreen attribute value - int type
*/
int Outlist::getToscreen() { return this->toscreen; }

//! get tofile attribute value
/*!
\return The tofile attribute value - int type
*/
int Outlist::getTofile() { return this->tofile; }

//! get thename attribute value
/*!
\return The thename attribute value - Array<char> * type
*/
Array<char> *Outlist::getThename() { return this->thename; }

// c++ attribute mutator methods
//! set toscreen attribute value
/*!
\param toscreenParam - the toscreen attribute value - int type
*/
void Outlist::setToscreen(int toscreenParam) { this->toscreen = toscreenParam; }

//! set tofile attribute value
/*!
\param tofileParam - the tofile attribute value - int type
*/
void Outlist::setTofile(int tofileParam) { this->tofile = tofileParam; }

//! set thename attribute value
/*!
\param thenameParam - the thename attribute value - Array<char> * type
*/
void Outlist::setThename(Array<char> *thenameParam) { this->thename = thenameParam; }

//######################################### array get/set #########################################

//! get thename attribute array element
/*!
\param i - indice - int type
\return The thename attribute value - char type
*/
char Outlist::getThenameElem(int i){
  int index[1] = {i};
  return this->thename->get(index);
}

//! set thename attribute array element
/*!
\param i - indice - int type
\param value - the value to be set - char type
*/
void Outlist::setThenameElem(int i, char value) {
  int index[1] = {i};
  this->thename->set(index, value);
  //std::cout<<i<<" "<<value<<std::endl; 
}

//#################################################################################################
//#################################################################################################

// Outlist instance pointer
Outlist *outlist;

// Fortran / C wrapper functions
extern "C" {
//! create Outlist instance
void createoutlist_() { outlist = new Outlist(); }

//! delete Outlist instance
void destroyoutlist_() { delete outlist; }

//! create thename array attribute
/*!
\param i - indice - int type
*/
void createoutlistthename_(int* i) { outlist->createThename(*i); }

//! delete thename array
void destroyoutlistthename_() { delete outlist->getThename(); }

//########################################## get/set ##############################################

//! get toscreen attribute value
/*!
\return The toscreen attribute value - int type
*/
int getoutlisttoscreen_() { return outlist->getToscreen(); }

//! set toscreen attribute value
/*!
\param toscreenParam - the toscreen attribute value - int type
*/
void setoutlisttoscreen_(int *toscreenParam) { outlist->setToscreen(*toscreenParam); }

//! get tofile attribute value
/*!
\return The tofile attribute value - int type
*/
int getoutlisttofile_() { return outlist->getTofile(); }

//! set tofile attribute value
/*!
\param tofileParam - the tofile attribute value - int type
*/
void setoutlisttofile_(int *tofileParam) { outlist->setTofile(*tofileParam); }

//##################################### array get/set #############################################

//! set thename array value
/*!
\param i - indice - int pointer type
\param value - value to set - char pointer type
*/
void setoutlistthename_(int *i, char *value) { outlist->setThenameElem(*i, *value); }

//! get thename array value
/*!
\param i indice - int pointer type
\param value - value to get - char pointer type
*/
void getoutlistthename_(int *i, char *value) { *value = outlist->getThenameElem(*i); }

}
