#ifndef OUTLIST_H
#define OUTLIST_H
#include "Array.hpp"
                     
class Outlist {

private:

  /*!
   * toscreen
   */
  int toscreen;

  /*!
   * tofile
   */
  int tofile;

  /*!
   * thename
   */
  Array<char> * thename;

public:

  //! default class constructor 
  Outlist();

  //! default class destructor 
  ~Outlist();

  // c++ array creation methods
  void createThename(int i);

  // c++ attribute accessor methods
  int getToscreen();
  int getTofile();
  Array<char> * getThename();

  // c++ attribute mutator methods
  void setToscreen(int toscreenParam);
  void setTofile(int tofileParam);
  void setThename(Array<char> * thenameParam);

  // c++ array accessor methods
  char getThenameElem(int i);

  // c++ array mutator methods
  void setThenameElem(int i, char value);

};
#endif
