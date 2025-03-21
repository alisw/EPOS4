#ifndef EOSP_H
#define EOSP_H
#include "Array.hpp"
                     
class Eosp {
private:
  /*!
   * oEeos
   */
  float oEeos;

public:
  //! default class constructor 
  Eosp();

  //! default class destructor 
  ~Eosp();

  // get value from fortran common block
  void getFortranOEeos();

  // set value from fortran common block
  void setFortranOEeos();

  // c++ attribute accessor methods
  float getOEeos();

  // c++ attribute mutator methods
  void setOEeos(float oEeosParam);

};
#endif
