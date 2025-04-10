#include "Event.hpp"

//! default class constructor 
Event::Event() {}

//! default class destructor 
Event::~Event() {
}

// map the fortran common block copcnt
extern "C" {
extern struct {
  int iopcnt;
  int jcentrality;
} copcnt_;
}

// get value from fortran common block
void Event::getFortranIopcnt() { this->iopcnt = copcnt_.iopcnt; }

// set value from fortran common block
void Event::setFortranIopcnt() { copcnt_.iopcnt = this->iopcnt; }

// c++ attribute accessor methods
//! get iopcnt attribute value
/*!
\return The iopcnt attribute value - int type
*/
int Event::getIopcnt() {
  // get value from fortran common block
  this->getFortranIopcnt();
  return this->iopcnt;
}

// c++ attribute mutator methods
//! set iopcnt attribute value
/*!
\param iopcntParam - the iopcnt attribute value - int type
*/
void Event::setIopcnt(int iopcntParam) {
  this->iopcnt = iopcntParam;
  // set value from fortran common block
  this->setFortranIopcnt();
}

// Event instance pointer
Event *event;

// Fortran / C wrapper functions
extern "C" {
//! create Event instance
void createevent_() { event = new Event(); }

//! delete Event instance
void destroyevent_() { delete event; }

//! get iopcnt attribute value
/*!
\return The iopcnt attribute value - int type
*/
int geteventiopcnt_() { return event->getIopcnt(); }

//! set iopcnt attribute value
/*!
\param iopcntParam - the iopcnt attribute value - int type
*/
void seteventiopcnt_(int *iopcntParam) { event->setIopcnt(*iopcntParam); }

}
