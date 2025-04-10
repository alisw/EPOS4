#ifndef CCOTI_H
#define CCOTI_H
#include "Array.hpp"
                     
class Ccoti {
private:
  /*!
   * computing time
   */
  Array<float> * coti;

public:
  //! default class constructor 
  Ccoti();

  //! default class destructor 
  ~Ccoti();

  // c++ attribute accessor methods
  Array<float> * getCoti();

  // c++ attribute mutator methods
  void setCoti(Array<float> * cotiParam);

  // c++ array creation methods
  void createCoti(int x, int y);

  // c++ array accessor methods
  float getCotiValue(int x, int y);

  // c++ array mutator methods
  void setCotiValue(int x, int y, float value);

  // c++ array initialize method
  void initializeCoti(float value);

};
#endif
