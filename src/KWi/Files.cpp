#include "Files.hpp"

//! default class constructor 
Files::Files() {}

//! default class destructor 
Files::~Files() {
}

// map the fortran common block nfname
extern "C" {
extern struct {
  int nfnch;
  int nfnhi;
  int nfndt;
  int nfnii;
  int nfnid;
  int nfnie;
  int nfnrj;
  int nfnmt;
  int nfngrv;
  int nfnnx;
  int nfncp;
  int nfncs;
  int nfndr;
  int nfnio;
} nfname_;
}

// get value from fortran common block
void Files::getFortranNfnio() { this->nfnio = nfname_.nfnio; }

// set value from fortran common block
void Files::setFortranNfnio() { nfname_.nfnio = this->nfnio; }

// map the fortran common block files
extern "C" {
extern struct {
  int ifop;
  int ifmt;
  int ifch;
  int ifcx;
  int ifhi;
  int ifdt;
  int ifcp;
  int ifdr;
  int ifio;
} files_;
}

// get value from fortran common block
void Files::getFortranIfmt() { this->ifmt = files_.ifmt; }

// set value from fortran common block
void Files::setFortranIfmt() { files_.ifmt = this->ifmt; }

// c++ attribute accessor methods
//! get nfnio attribute value
/*!
\return The nfnio attribute value - int type
*/
int Files::getNfnio() {
  // get value from fortran common block
  this->getFortranNfnio();
  return this->nfnio;
}

//! get ifmt attribute value
/*!
\return The ifmt attribute value - int type
*/
int Files::getIfmt() {
  // get value from fortran common block
  this->getFortranIfmt();
  return this->ifmt;
}

//! get fnio attribute value
/*!
\return The fnio attribute value - Array<char> * type
*/
Array<char> *Files::getFnio() { return this->fnio; }

// c++ attribute mutator methods
//! set nfnio attribute value
/*!
\param nfnioParam - the nfnio attribute value - int type
*/
void Files::setNfnio(int nfnioParam) {
  this->nfnio = nfnioParam;
  // set value from fortran common block
  this->setFortranNfnio();
}

//! set ifmt attribute value
/*!
\param ifmtParam - the ifmt attribute value - int type
*/
void Files::setIfmt(int ifmtParam) {
  this->ifmt = ifmtParam;
  // set value from fortran common block
  this->setFortranIfmt();
}

//! set fnio attribute value
/*!
\param fnioParam - the fnio attribute value - Array<char> * type
*/
void Files::setFnio(Array<char> *fnioParam) { this->fnio = fnioParam; }

// c++ array attribute methods
//! create fnio attribute array of char type
/*!
\param i - indice - int type
*/
void Files::createFnio(int i) {
  int dim = 1;
  int shape[1] = {i};
  this->fnio = new Array<char>(dim, shape);
}

//! get fnio attribute array value
/*!
\param i - indice - int type
\return The fnio attribute value - char type
*/
char Files::getFnioValue(int i){
  int index[1] = {i};
  return this->fnio->get(index);
}

//! set fnio attribute array value
/*!
\param i - indice - int type
\param value - the value to be set - char type
*/
void Files::setFnioValue(int i, char value) {
  int index[1] = {i};
  this->fnio->set(index, value);
}

//! create fnio initialize method
void Files::initializeFnio(char value) {
  this->fnio->initialize(value);
}

// Files instance pointer
Files *files;

// Fortran / C wrapper functions
extern "C" {
//! create Files instance
void createfiles_() { files = new Files(); }

//! delete Files instance
void destroyfiles_() { delete files; }

//! get nfnio attribute value
/*!
\return The nfnio attribute value - int type
*/
int getfilesnfnio_() { return files->getNfnio(); }

//! set nfnio attribute value
/*!
\param nfnioParam - the nfnio attribute value - int type
*/
void setfilesnfnio_(int *nfnioParam) { files->setNfnio(*nfnioParam); }

//! get ifmt attribute value
/*!
\return The ifmt attribute value - int type
*/
int getfilesifmt_() { return files->getIfmt(); }

//! set ifmt attribute value
/*!
\param ifmtParam - the ifmt attribute value - int type
*/
void setfilesifmt_(int *ifmtParam) { files->setIfmt(*ifmtParam); }

//! create fnio array attribute
/*!\param i - indice - int type
*/void createfilesfnio_(int* i) {
  files->createFnio(*i);
}

//! initialize fnio array attribute
void initializefilesfnio_(char *value) {
  files->initializeFnio(*value);
}

//! set fnio array value
/*!
\param i - indice - int pointer type
\param value - value to set - char pointer type
*/
void setfilesfnio_(int *i, char *value) {
  files->setFnioValue(*i, *value);
}

//! get fnio array value
/*!
\param i indice - int pointer type
\param value - value to get - char pointer type
*/
void getfilesfnio_(int *i, char *value) {
  *value = files->getFnioValue(*i);
}

//! delete fnio array
void destroyfilesfnio_() { delete files->getFnio(); }

}
