#ifndef ACCUM_H
#define ACCUM_H
#include "Array.hpp"
                     
class Accum {
private:
  /*!
   * imsg
   */
  int imsg;
  /*!
   * computing time
   */
  Array<int> * jerr;
  /*!
   * ntevt
   */
  int ntevt;
  /*!
   * nrevt
   */
  int nrevt;

public:
  //! default class constructor 
  Accum();

  //! default class destructor 
  ~Accum();

  // get value from fortran common block
  void getFortranImsg();
  void getFortranNtevt();
  void getFortranNrevt();

  // set value from fortran common block
  void setFortranImsg();
  void setFortranNtevt();
  void setFortranNrevt();

  // c++ attribute accessor methods
  int getImsg();
  Array<int> * getJerr();
  int getNtevt();
  int getNrevt();

  // c++ attribute mutator methods
  void setImsg(int imsgParam);
  void setJerr(Array<int> * jerrParam);
  void setNtevt(int ntevtParam);
  void setNrevt(int nrevtParam);

  // c++ array creation methods
  void createJerr(int i);

  // c++ array accessor methods
  int getJerrValue(int i);

  // c++ array mutator methods
  void setJerrValue(int i, int value);

  // c++ array initialize method
  void initializeJerr(int value);

};
#endif
