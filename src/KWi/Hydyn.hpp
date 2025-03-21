#ifndef HYDYN_H
#define HYDYN_H
#include "Array.hpp"
                     
class Hydyn {
private:
  /*!
   * nxhy
   */
  int nxhy;
  /*!
   * nyhy
   */
  int nyhy;
  /*!
   * nzhy
   */
  int nzhy;
  /*!
   * ntauhy
   */
  int ntauhy;
  /*!
   * tau step for hydro
   */
  float dtauhy;
  /*!
   * x min for hydro output
   */
  float xminhy;
  /*!
   * x max for hydro output
   */
  float xmaxhy;
  /*!
   * y min for hydro output
   */
  float yminhy;
  /*!
   * y max for hydro output
   */
  float ymaxhy;
  /*!
   * eta min for hydro output
   */
  float zminhy;
  /*!
   * eta max for hydro output
   */
  float zmaxhy;
  /*!
   * tfrout
   */
  float tfrout;
  /*!
   * ioeos
   */
  int ioeos;
  /*!
   * iozerof
   */
  int iozerof;
  /*!
   * ihlle
   */
  int ihlle;
  /*!
   * jhlle
   */
  int jhlle;
  /*!
   * ntaumx
   */
  int ntaumx;
  /*!
   * tfo
   */
  float tfo;
  /*!
   * activate envar (1) or not (0) !1 = very slow! only test
   */
  int ienvar;
  /*!
   * ifaahlle
   */
  int ifaahlle;
  /*!
   * ifathlle
   */
  int ifathlle;
  /*!
   * ifazhlle
   */
  int ifazhlle;
  /*!
   * epsfin
   */
  float epsfin;
  /*!
   * tauhy
   */
  Array<float> * tauhy;
  /*!
   * eetau
   */
  float eetau;
  /*!
   * Energy from TransferHlle (sparce grid)
   */
  float eetau2;
  /*!
   * emx
   */
  float emx;
  /*!
   * xmx
   */
  float xmx;
  /*!
   * ymx
   */
  float ymx;
  /*!
   * zmx
   */
  float zmx;
  /*!
   * mmx
   */
  int mmx;
  /*!
   * mmy
   */
  int mmy;
  /*!
   * mmz
   */
  int mmz;
  /*!
   * ratioeex
   */
  Array<float> * ratioeex;
  /*!
   * velc
   */
  Array<double> * velc;
  /*!
   * epsc
   */
  Array<double> * epsc;
  /*!
   * sigc
   */
  Array<double> * sigc;
  /*!
   * barc
   */
  Array<double> * barc;

public:
  //! default class constructor 
  Hydyn();

  Hydyn (int ntau);

  //! default class destructor 
  ~Hydyn();

  // get value from fortran common block
  void getFortranNxhy();
  void getFortranNyhy();
  void getFortranNzhy();
  void getFortranNtauhy();

  // set value from fortran common block
  void setFortranNxhy();
  void setFortranNyhy();
  void setFortranNzhy();
  void setFortranNtauhy();

  // get value from fortran common block
  void getFortranDtauhy();

  // set value from fortran common block
  void setFortranDtauhy();

  // get value from fortran common block
  void getFortranXminhy();
  void getFortranXmaxhy();
  void getFortranYminhy();
  void getFortranYmaxhy();
  void getFortranZminhy();
  void getFortranZmaxhy();

  // set value from fortran common block
  void setFortranXminhy();
  void setFortranXmaxhy();
  void setFortranYminhy();
  void setFortranYmaxhy();
  void setFortranZminhy();
  void setFortranZmaxhy();

  // get value from fortran common block
  void getFortranTfrout();

  // set value from fortran common block
  void setFortranTfrout();

  // get value from fortran common block
  void getFortranIoeos();
  void getFortranIozerof();

  // set value from fortran common block
  void setFortranIoeos();
  void setFortranIozerof();

  // get value from fortran common block
  void getFortranIhlle();
  void getFortranJhlle();
  void getFortranNtaumx();

  // set value from fortran common block
  void setFortranIhlle();
  void setFortranJhlle();
  void setFortranNtaumx();

  // get value from fortran common block
  void getFortranTfo();

  // set value from fortran common block
  void setFortranTfo();

  // get value from fortran common block
  void getFortranIenvar();

  // set value from fortran common block
  void setFortranIenvar();

  // get value from fortran common block
  void getFortranIfaahlle();
  void getFortranIfathlle();
  void getFortranIfazhlle();

  // set value from fortran common block
  void setFortranIfaahlle();
  void setFortranIfathlle();
  void setFortranIfazhlle();

  // get value from fortran common block
  void getFortranEpsfin();

  // set value from fortran common block
  void setFortranEpsfin();

  // c++ attribute accessor methods
  int getNxhy();
  int getNyhy();
  int getNzhy();
  int getNtauhy();
  float getDtauhy();
  float getXminhy();
  float getXmaxhy();
  float getYminhy();
  float getYmaxhy();
  float getZminhy();
  float getZmaxhy();
  float getTfrout();
  int getIoeos();
  int getIozerof();
  int getIhlle();
  int getJhlle();
  int getNtaumx();
  float getTfo();
  int getIenvar();
  int getIfaahlle();
  int getIfathlle();
  int getIfazhlle();
  float getEpsfin();
  Array<float> * getTauhy();
  float getEetau();
  float getEetau2();
  float getEmx();
  float getXmx();
  float getYmx();
  float getZmx();
  int getMmx();
  int getMmy();
  int getMmz();
  Array<float> * getRatioeex();
  Array<double> * getVelc();
  Array<double> * getEpsc();
  Array<double> * getSigc();
  Array<double> * getBarc();

  // c++ attribute mutator methods
  void setNxhy(int nxhyParam);
  void setNyhy(int nyhyParam);
  void setNzhy(int nzhyParam);
  void setNtauhy(int ntauhyParam);
  void setDtauhy(float dtauhyParam);
  void setXminhy(float xminhyParam);
  void setXmaxhy(float xmaxhyParam);
  void setYminhy(float yminhyParam);
  void setYmaxhy(float ymaxhyParam);
  void setZminhy(float zminhyParam);
  void setZmaxhy(float zmaxhyParam);
  void setTfrout(float tfroutParam);
  void setIoeos(int ioeosParam);
  void setIozerof(int iozerofParam);
  void setIhlle(int ihlleParam);
  void setJhlle(int jhlleParam);
  void setNtaumx(int ntaumxParam);
  void setTfo(float tfoParam);
  void setIenvar(int ienvarParam);
  void setIfaahlle(int ifaahlleParam);
  void setIfathlle(int ifathlleParam);
  void setIfazhlle(int ifazhlleParam);
  void setEpsfin(float epsfinParam);
  void setTauhy(Array<float> * tauhyParam);
  void setEetau(float eetauParam);
  void setEetau2(float eetau2Param);
  void setEmx(float emxParam);
  void setXmx(float xmxParam);
  void setYmx(float ymxParam);
  void setZmx(float zmxParam);
  void setMmx(int mmxParam);
  void setMmy(int mmyParam);
  void setMmz(int mmzParam);
  void setRatioeex(Array<float> * ratioeexParam);
  void setVelc(Array<double> * velcParam);
  void setEpsc(Array<double> * epscParam);
  void setSigc(Array<double> * sigcParam);
  void setBarc(Array<double> * barcParam);

  // c++ multiple attribute mutator methods
  void setDim(int nxhyParam, int nyhyParam, int nzhyParam, int ntauhyParam);
  void setBound(float xminhyParam, float xmaxhyParam, float yminhyParam, float ymaxhyParam, float zminhyParam, float zmaxhyParam);
  void setIo(int ioeosParam, int iozerofParam);

  // c++ multiple attribute accessor methods
  void getDim(int* nxhyParam, int* nyhyParam, int* nzhyParam, int* ntauhyParam);
  void getBound(float* xminhyParam, float* xmaxhyParam, float* yminhyParam, float* ymaxhyParam, float* zminhyParam, float* zmaxhyParam);
  void getIo(int* ioeosParam, int* iozerofParam);

  // c++ array creation methods
  void createTauhy(int ntau);
  void createRatioeex(int ntau);
  void createVelc(int n, int neta, int ntau, int nx, int ny);
  void createEpsc(int neta, int ntau, int nx, int ny);
  void createSigc(int neta, int ntau, int nx, int ny);
  void createBarc(int n, int neta, int ntau, int nx, int ny);

  // c++ array accessor methods
  float getTauhyValue(int ntau);
  float getRatioeexValue(int ntau);
  double getVelcValue(int n, int neta, int ntau, int nx, int ny);
  double getEpscValue(int neta, int ntau, int nx, int ny);
  double getSigcValue(int neta, int ntau, int nx, int ny);
  double getBarcValue(int n, int neta, int ntau, int nx, int ny);

  // c++ array mutator methods
  void setTauhyValue(int ntau, float value);
  void setRatioeexValue(int ntau, float value);
  void setVelcValue(int n, int neta, int ntau, int nx, int ny, double value);
  void setEpscValue(int neta, int ntau, int nx, int ny, double value);
  void setSigcValue(int neta, int ntau, int nx, int ny, double value);
  void setBarcValue(int n, int neta, int ntau, int nx, int ny, double value);

  // c++ array initialize method
  void initializeTauhy(const float value);
  void initializeRatioeex(const float value);
  void initializeVelc(const double value);
  void initializeEpsc(const double value);
  void initializeSigc(const double value);
  void initializeBarc(const double value);

};
#endif
