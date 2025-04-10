#ifndef CORE_H
#define CORE_H
#include "Array.hpp"
                     
class Core {
private:
  /*!
   * read table for initial condition stuff
   */
  int icotabr;

public:
  //! default class constructor 
  Core();

  //! default class destructor 
  ~Core();

  // get value from fortran common block
  void getFortranIcotabr();

  // set value from fortran common block
  void setFortranIcotabr();

  // c++ attribute accessor methods
  int getIcotabr();

  // c++ attribute mutator methods
  void setIcotabr(int icotabrParam);

};
#endif
