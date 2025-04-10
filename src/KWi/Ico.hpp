#ifndef ICO_H
#define ICO_H
#include "Array.hpp"

/* --------------------------------------------------------------------
  initial conditions for hydro
---------------------------------------------------------------------*/

/* --------------------------------------------------------------------
  We think in terms of bins:
  nxico,nyico,nzico: number of bins
  xminico, yminico,zminico: lower end of 1st bin
  xmaxico,ymaxico,zmaxico: upper end of last bin
  Concerning the mean value for the bins numbers i,j,k:
  x_i=xminico+(i-0.5)*(xmaxico-xminico)/nxico
  y_j=yminico+(j-0.5)*(ymaxico-yminico)/nyico
  z_k=zminico+(k-0.5)*(zmaxico-zminico)/nzico
                                                                                                                                                              
  ATTENTION:        z means eta_s !!!!!!
---------------------------------------------------------------------*/

/*---------------------------------------------------------------------
                                                                                                                                    
  Primary tables (from string segments):
----------------------------------------
  /Ico1/  IcoT(ip1,ip2,ix,iy,iz) energy-momentum Tensor
  /Ico2/  IcoC(if,ip,ix,iy,iz)   flavor current 4-vector
                                                                                                                                                   
  Secondary tables (after diadonalizing T):
-------------------------------------------
  /Ico3/  IcoE(ix,iy,iz)         energy-density
  /Ico4/  IcoV(1,ix,iy,iz)       x-velocity
          IcoV(2,ix,iy,iz)       y-velocity
          IcoV(3,ix,iy,iz)       eta-velocity (z-velocity/tau)
  /Ico5/  IcoF(1,ix,iy,iz)       net flavor density of up
          IcoF(2,ix,iy,iz)       net flavor density of down
          IcoF(3,ix,iy,iz)       net flavor density of strange
----------------------------------------------------------------------
  ip, ip1,ip2 .... dirac indices
  ix,iy,iz ....... transverse coordinates + pseudo-rapidity eta
  if ............. flavor
----------------------------------------------------------------------
  To read in secondary tables, see subroutine IniCon:
  read(97,*) ...
---------------------------------------------------------------------*/


class Ico {
private:
  /*!
   * xmin for initial condition calculation
   */
  float xmin;
  /*!
   * xmax for initial condition calculation
   */
  float xmax;
  /*!
   * ymin for initial condition calculation
   */
  float ymin;
  /*!
   * ymax for initial condition calculation
   */
  float ymax;
  /*!
   * zmin for initial condition calculation
   */
  float zmin;
  /*!
   * zmax for initial condition calculation
   */
  float zmax;
  /*!
   * nx
   */
  int nx;
  /*!
   * ny
   */
  int ny;
  /*!
   * nz
   */
  int nz;
  /*!
   * E from Ico
   */
  float ee1ico;
  /*!
   * eistico
   */
  float eistico;
  /*!
   * ee1hll
   */
  float ee1hll;
  /*!
   * ichkengy
   */
  int ichkengy;
  /*!
   * esollxx
   */
  float esollxx;
  /*!
   * eistxx 
   */
  float eistxx;
  /*!
   * energy density
   */
  Array<float> * energyDensity;
  /*!
   * flavor density
   */
  Array<float> * flavorDensity;
  /*!
   * velocity
   */
  Array<float> * velocity;
  /*!
   * energy momentum tensor
   */
  Array<double> * energyMomentumTensor;
  /*!
   * flavor current vector
   */
  Array<double> * flavorCurrentVector;

public:
  //! default class constructor 
  Ico();

  //! default class destructor 
  ~Ico();

  // get value from fortran common block
  void getFortranEsollxx();
  void getFortranEistxx();

  // set value from fortran common block
  void setFortranEsollxx();
  void setFortranEistxx();

  // c++ attribute accessor methods
  float getXmin();
  float getXmax();
  float getYmin();
  float getYmax();
  float getZmin();
  float getZmax();
  int getNx();
  int getNy();
  int getNz();
  float getEe1ico();
  float getEistico();
  float getEe1hll();
  int getIchkengy();
  float getEsollxx();
  float getEistxx();
  Array<float> * getEnergyDensity();
  Array<float> * getFlavorDensity();
  Array<float> * getVelocity();
  Array<double> * getEnergyMomentumTensor();
  Array<double> * getFlavorCurrentVector();

  // c++ attribute mutator methods
  void setXmin(float xminParam);
  void setXmax(float xmaxParam);
  void setYmin(float yminParam);
  void setYmax(float ymaxParam);
  void setZmin(float zminParam);
  void setZmax(float zmaxParam);
  void setNx(int nxParam);
  void setNy(int nyParam);
  void setNz(int nzParam);
  void setEe1ico(float ee1icoParam);
  void setEistico(float eisticoParam);
  void setEe1hll(float ee1hllParam);
  void setIchkengy(int ichkengyParam);
  void setEsollxx(float esollxxParam);
  void setEistxx(float eistxxParam);
  void setEnergyDensity(Array<float> * energyDensityParam);
  void setFlavorDensity(Array<float> * flavorDensityParam);
  void setVelocity(Array<float> * velocityParam);
  void setEnergyMomentumTensor(Array<double> * energyMomentumTensorParam);
  void setFlavorCurrentVector(Array<double> * flavorCurrentVectorParam);

  // c++ multiple attribute mutator methods
  void setDim(int nxParam, int nyParam, int nzParam);
  void setBound(float xminParam, float xmaxParam, float yminParam, float ymaxParam, float zminParam, float zmaxParam);

  // c++ multiple attribute accessor methods
  void getDim(int* nxParam, int* nyParam, int* nzParam);
  void getBound(float* xminParam, float* xmaxParam, float* yminParam, float* ymaxParam, float* zminParam, float* zmaxParam);

  // c++ array creation methods
  void createEnergyDensity(int x, int y, int z);
  void createFlavorDensity(int n, int x, int y, int z);
  void createVelocity(int n, int x, int y, int z);
  void createEnergyMomentumTensor(int ip1, int ip2, int x, int y, int z);
  void createFlavorCurrentVector(int f, int ip, int x, int y, int z);

  // c++ array accessor methods
  float getEnergyDensityValue(int x, int y, int z);
  float getFlavorDensityValue(int n, int x, int y, int z);
  float getVelocityValue(int n, int x, int y, int z);
  double getEnergyMomentumTensorValue(int ip1, int ip2, int x, int y, int z);
  double getFlavorCurrentVectorValue(int f, int ip, int x, int y, int z);

  // c++ array mutator methods
  void setEnergyDensityValue(int x, int y, int z, float value);
  void setFlavorDensityValue(int n, int x, int y, int z, float value);
  void setVelocityValue(int n, int x, int y, int z, float value);
  void setEnergyMomentumTensorValue(int ip1, int ip2, int x, int y, int z, double value);
  void setFlavorCurrentVectorValue(int f, int ip, int x, int y, int z, double value);

  // c++ array initialize method
  void initializeEnergyDensity(float value);
  void initializeFlavorDensity(float value);
  void initializeVelocity(float value);
  void initializeEnergyMomentumTensor(double value);
  void initializeFlavorCurrentVector(double value);

};
#endif
