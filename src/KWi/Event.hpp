#ifndef EVENT_H
#define EVENT_H
#include "Array.hpp"
                     
class Event {
private:
  /*!
   * iopcnt
   */
  int iopcnt;

public:
  //! default class constructor 
  Event();

  //! default class destructor 
  ~Event();

  // get value from fortran common block
  void getFortranIopcnt();

  // set value from fortran common block
  void setFortranIopcnt();

  // c++ attribute accessor methods
  int getIopcnt();

  // c++ attribute mutator methods
  void setIopcnt(int iopcntParam);

};
#endif
