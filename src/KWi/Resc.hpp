#ifndef RESC_H
#define RESC_H
#include "Array.hpp"
                     
class Resc {
private:
  /*!
   * tau for core-corona procedure
   */
  float tauZero;
  /*!
   * hydro tau step / 2 for tau >  tauone
   */
  float tauOne;
  /*!
   * hydro tau step / 4 for tau >  tautwo
   */
  float tauTwo;
  /*!
   * tau three
   */
  float tauThree;
  /*!
   * upper limit for tau - tauzer for hydro
   */
  float tauUp;
  /*!
   * etaos
   */
  float etaos;
  /*!
   * zetaos
   */
  float zetaos;

public:
  //! default class constructor 
  Resc();

  //! default class destructor 
  ~Resc();

  // get value from fortran common block
  void getFortranTauZero();
  void getFortranTauOne();
  void getFortranTauTwo();
  void getFortranTauUp();
  void getFortranTauThree();

  // set value from fortran common block
  void setFortranTauZero();
  void setFortranTauOne();
  void setFortranTauTwo();
  void setFortranTauUp();
  void setFortranTauThree();

  // get value from fortran common block
  void getFortranEtaos();
  void getFortranZetaos();

  // set value from fortran common block
  void setFortranEtaos();
  void setFortranZetaos();

  // c++ attribute accessor methods
  float getTauZero();
  float getTauOne();
  float getTauTwo();
  float getTauThree();
  float getTauUp();
  float getEtaos();
  float getZetaos();

  // c++ attribute mutator methods
  void setTauZero(float tauZeroParam);
  void setTauOne(float tauOneParam);
  void setTauTwo(float tauTwoParam);
  void setTauThree(float tauThreeParam);
  void setTauUp(float tauUpParam);
  void setEtaos(float etaosParam);
  void setZetaos(float zetaosParam);

  // c++ multiple attribute mutator methods
  void setTau(float tauZeroParam, float tauOneParam, float tauTwoParam, float tauThreeParam, float tauUpParam);

  // c++ multiple attribute accessor methods
  void getTau(float* tauZeroParam, float* tauOneParam, float* tauTwoParam, float* tauThreeParam, float* tauUpParam);

};
#endif
