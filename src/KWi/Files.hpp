#ifndef FILES_H
#define FILES_H
#include "Array.hpp"
                     
class Files {
private:
  /*!
   * nfnio
   */
  int nfnio;
  /*!
   * std output unit
   */
  int ifmt;
  /*!
   * fnio
   */
  Array<char> * fnio;

public:
  //! default class constructor 
  Files();

  //! default class destructor 
  ~Files();

  // get value from fortran common block
  void getFortranNfnio();

  // set value from fortran common block
  void setFortranNfnio();

  // get value from fortran common block
  void getFortranIfmt();

  // set value from fortran common block
  void setFortranIfmt();

  // c++ attribute accessor methods
  int getNfnio();
  int getIfmt();
  Array<char> * getFnio();

  // c++ attribute mutator methods
  void setNfnio(int nfnioParam);
  void setIfmt(int ifmtParam);
  void setFnio(Array<char> * fnioParam);

  // c++ array creation methods
  void createFnio(int i);

  // c++ array accessor methods
  char getFnioValue(int i);

  // c++ array mutator methods
  void setFnioValue(int i, char value);

  // c++ array initialize method
  void initializeFnio(char value);

};
#endif
