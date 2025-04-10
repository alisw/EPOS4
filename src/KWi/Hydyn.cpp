#include "Hydyn.hpp"

//! default class constructor 
Hydyn::Hydyn() {}

//! default class destructor 
Hydyn::~Hydyn() {
}

// map the fortran common block hocoxyeta
extern "C" {
extern struct {
  int nxhy;
  int nyhy;
  int nzhy;
  int ntauhy;
} hocoxyeta_;
}

// get value from fortran common block
void Hydyn::getFortranNxhy() { this->nxhy = hocoxyeta_.nxhy; }

void Hydyn::getFortranNyhy() { this->nyhy = hocoxyeta_.nyhy; }

void Hydyn::getFortranNzhy() { this->nzhy = hocoxyeta_.nzhy; }

void Hydyn::getFortranNtauhy() { this->ntauhy = hocoxyeta_.ntauhy; }

// set value from fortran common block
void Hydyn::setFortranNxhy() { hocoxyeta_.nxhy = this->nxhy; }

void Hydyn::setFortranNyhy() { hocoxyeta_.nyhy = this->nyhy; }

void Hydyn::setFortranNzhy() { hocoxyeta_.nzhy = this->nzhy; }

void Hydyn::setFortranNtauhy() { hocoxyeta_.ntauhy = this->ntauhy; }

// map the fortran common block hocotau2
extern "C" {
extern struct {
  float dtauhy;
} hocotau2_;
}

// get value from fortran common block
void Hydyn::getFortranDtauhy() { this->dtauhy = hocotau2_.dtauhy; }

// set value from fortran common block
void Hydyn::setFortranDtauhy() { hocotau2_.dtauhy = this->dtauhy; }

// map the fortran common block hocobound
extern "C" {
extern struct {
  float xminhy;
  float xmaxhy;
  float yminhy;
  float ymaxhy;
  float zminhy;
  float zmaxhy;
} hocobound_;
}

// get value from fortran common block
void Hydyn::getFortranXminhy() { this->xminhy = hocobound_.xminhy; }

void Hydyn::getFortranXmaxhy() { this->xmaxhy = hocobound_.xmaxhy; }

void Hydyn::getFortranYminhy() { this->yminhy = hocobound_.yminhy; }

void Hydyn::getFortranYmaxhy() { this->ymaxhy = hocobound_.ymaxhy; }

void Hydyn::getFortranZminhy() { this->zminhy = hocobound_.zminhy; }

void Hydyn::getFortranZmaxhy() { this->zmaxhy = hocobound_.zmaxhy; }

// set value from fortran common block
void Hydyn::setFortranXminhy() { hocobound_.xminhy = this->xminhy; }

void Hydyn::setFortranXmaxhy() { hocobound_.xmaxhy = this->xmaxhy; }

void Hydyn::setFortranYminhy() { hocobound_.yminhy = this->yminhy; }

void Hydyn::setFortranYmaxhy() { hocobound_.ymaxhy = this->ymaxhy; }

void Hydyn::setFortranZminhy() { hocobound_.zminhy = this->zminhy; }

void Hydyn::setFortranZmaxhy() { hocobound_.zmaxhy = this->zmaxhy; }

// map the fortran common block hoco13
extern "C" {
extern struct {
  float tfrout;
} hoco13_;
}

// get value from fortran common block
void Hydyn::getFortranTfrout() { this->tfrout = hoco13_.tfrout; }

// set value from fortran common block
void Hydyn::setFortranTfrout() { hoco13_.tfrout = this->tfrout; }

// map the fortran common block hydr3
extern "C" {
extern struct {
  int ioeos;
  int iozerof;
} hydr3_;
}

// get value from fortran common block
void Hydyn::getFortranIoeos() { this->ioeos = hydr3_.ioeos; }

void Hydyn::getFortranIozerof() { this->iozerof = hydr3_.iozerof; }

// set value from fortran common block
void Hydyn::setFortranIoeos() { hydr3_.ioeos = this->ioeos; }

void Hydyn::setFortranIozerof() { hydr3_.iozerof = this->iozerof; }

// map the fortran common block hlle1
extern "C" {
extern struct {
  int ihlle;
  int jhlle;
  int ntaumx;
} hlle1_;
}

// get value from fortran common block
void Hydyn::getFortranIhlle() { this->ihlle = hlle1_.ihlle; }

void Hydyn::getFortranJhlle() { this->jhlle = hlle1_.jhlle; }

void Hydyn::getFortranNtaumx() { this->ntaumx = hlle1_.ntaumx; }

// set value from fortran common block
void Hydyn::setFortranIhlle() { hlle1_.ihlle = this->ihlle; }

void Hydyn::setFortranJhlle() { hlle1_.jhlle = this->jhlle; }

void Hydyn::setFortranNtaumx() { hlle1_.ntaumx = this->ntaumx; }

// map the fortran common block ctfo
extern "C" {
extern struct {
  float tfo;
} ctfo_;
}

// get value from fortran common block
void Hydyn::getFortranTfo() { this->tfo = ctfo_.tfo; }

// set value from fortran common block
void Hydyn::setFortranTfo() { ctfo_.tfo = this->tfo; }

// map the fortran common block copt
extern "C" {
extern struct {
  int istat;
  int ienvar;
} copt_;
}

// get value from fortran common block
void Hydyn::getFortranIenvar() { this->ienvar = copt_.ienvar; }

// set value from fortran common block
void Hydyn::setFortranIenvar() { copt_.ienvar = this->ienvar; }

// map the fortran common block cifahlle
extern "C" {
extern struct {
  int ifaahlle;
  int ifathlle;
  int ifazhlle;
} cifahlle_;
}

// get value from fortran common block
void Hydyn::getFortranIfaahlle() { this->ifaahlle = cifahlle_.ifaahlle; }

void Hydyn::getFortranIfathlle() { this->ifathlle = cifahlle_.ifathlle; }

void Hydyn::getFortranIfazhlle() { this->ifazhlle = cifahlle_.ifazhlle; }

// set value from fortran common block
void Hydyn::setFortranIfaahlle() { cifahlle_.ifaahlle = this->ifaahlle; }

void Hydyn::setFortranIfathlle() { cifahlle_.ifathlle = this->ifathlle; }

void Hydyn::setFortranIfazhlle() { cifahlle_.ifazhlle = this->ifazhlle; }

// map the fortran common block cepsfin
extern "C" {
extern struct {
  float epsfin;
} cepsfin_;
}

// get value from fortran common block
void Hydyn::getFortranEpsfin() { this->epsfin = cepsfin_.epsfin; }

// set value from fortran common block
void Hydyn::setFortranEpsfin() { cepsfin_.epsfin = this->epsfin; }

// c++ attribute accessor methods
//! get nxhy attribute value
/*!
\return The nxhy attribute value - int type
*/
int Hydyn::getNxhy() {
  // get value from fortran common block
  this->getFortranNxhy();
  return this->nxhy;
}

//! get nyhy attribute value
/*!
\return The nyhy attribute value - int type
*/
int Hydyn::getNyhy() {
  // get value from fortran common block
  this->getFortranNyhy();
  return this->nyhy;
}

//! get nzhy attribute value
/*!
\return The nzhy attribute value - int type
*/
int Hydyn::getNzhy() {
  // get value from fortran common block
  this->getFortranNzhy();
  return this->nzhy;
}

//! get ntauhy attribute value
/*!
\return The ntauhy attribute value - int type
*/
int Hydyn::getNtauhy() {
  // get value from fortran common block
  this->getFortranNtauhy();
  return this->ntauhy;
}

//! get dtauhy attribute value
/*!
\return The dtauhy attribute value - float type
*/
float Hydyn::getDtauhy() {
  // get value from fortran common block
  this->getFortranDtauhy();
  return this->dtauhy;
}

//! get xminhy attribute value
/*!
\return The xminhy attribute value - float type
*/
float Hydyn::getXminhy() {
  // get value from fortran common block
  this->getFortranXminhy();
  return this->xminhy;
}

//! get xmaxhy attribute value
/*!
\return The xmaxhy attribute value - float type
*/
float Hydyn::getXmaxhy() {
  // get value from fortran common block
  this->getFortranXmaxhy();
  return this->xmaxhy;
}

//! get yminhy attribute value
/*!
\return The yminhy attribute value - float type
*/
float Hydyn::getYminhy() {
  // get value from fortran common block
  this->getFortranYminhy();
  return this->yminhy;
}

//! get ymaxhy attribute value
/*!
\return The ymaxhy attribute value - float type
*/
float Hydyn::getYmaxhy() {
  // get value from fortran common block
  this->getFortranYmaxhy();
  return this->ymaxhy;
}

//! get zminhy attribute value
/*!
\return The zminhy attribute value - float type
*/
float Hydyn::getZminhy() {
  // get value from fortran common block
  this->getFortranZminhy();
  return this->zminhy;
}

//! get zmaxhy attribute value
/*!
\return The zmaxhy attribute value - float type
*/
float Hydyn::getZmaxhy() {
  // get value from fortran common block
  this->getFortranZmaxhy();
  return this->zmaxhy;
}

//! get tfrout attribute value
/*!
\return The tfrout attribute value - float type
*/
float Hydyn::getTfrout() {
  // get value from fortran common block
  this->getFortranTfrout();
  return this->tfrout;
}

//! get ioeos attribute value
/*!
\return The ioeos attribute value - int type
*/
int Hydyn::getIoeos() {
  // get value from fortran common block
  this->getFortranIoeos();
  return this->ioeos;
}

//! get iozerof attribute value
/*!
\return The iozerof attribute value - int type
*/
int Hydyn::getIozerof() {
  // get value from fortran common block
  this->getFortranIozerof();
  return this->iozerof;
}

//! get ihlle attribute value
/*!
\return The ihlle attribute value - int type
*/
int Hydyn::getIhlle() {
  // get value from fortran common block
  this->getFortranIhlle();
  return this->ihlle;
}

//! get jhlle attribute value
/*!
\return The jhlle attribute value - int type
*/
int Hydyn::getJhlle() {
  // get value from fortran common block
  this->getFortranJhlle();
  return this->jhlle;
}

//! get ntaumx attribute value
/*!
\return The ntaumx attribute value - int type
*/
int Hydyn::getNtaumx() {
  // get value from fortran common block
  this->getFortranNtaumx();
  return this->ntaumx;
}

//! get tfo attribute value
/*!
\return The tfo attribute value - float type
*/
float Hydyn::getTfo() {
  // get value from fortran common block
  this->getFortranTfo();
  return this->tfo;
}

//! get ienvar attribute value
/*!
\return The ienvar attribute value - int type
*/
int Hydyn::getIenvar() {
  // get value from fortran common block
  this->getFortranIenvar();
  return this->ienvar;
}

//! get ifaahlle attribute value
/*!
\return The ifaahlle attribute value - int type
*/
int Hydyn::getIfaahlle() {
  // get value from fortran common block
  this->getFortranIfaahlle();
  return this->ifaahlle;
}

//! get ifathlle attribute value
/*!
\return The ifathlle attribute value - int type
*/
int Hydyn::getIfathlle() {
  // get value from fortran common block
  this->getFortranIfathlle();
  return this->ifathlle;
}

//! get ifazhlle attribute value
/*!
\return The ifazhlle attribute value - int type
*/
int Hydyn::getIfazhlle() {
  // get value from fortran common block
  this->getFortranIfazhlle();
  return this->ifazhlle;
}

//! get epsfin attribute value
/*!
\return The epsfin attribute value - float type
*/
float Hydyn::getEpsfin() {
  // get value from fortran common block
  this->getFortranEpsfin();
  return this->epsfin;
}

//! get tauhy attribute value
/*!
\return The tauhy attribute value - Array<float> * type
*/
Array<float> *Hydyn::getTauhy() { return this->tauhy; }

//! get eetau attribute value
/*!
\return The eetau attribute value - float type
*/
float Hydyn::getEetau() { return this->eetau; }

//! get eetau2 attribute value
/*!
\return The eetau2 attribute value - float type
*/
float Hydyn::getEetau2() { return this->eetau2; }

//! get emx attribute value
/*!
\return The emx attribute value - float type
*/
float Hydyn::getEmx() { return this->emx; }

//! get xmx attribute value
/*!
\return The xmx attribute value - float type
*/
float Hydyn::getXmx() { return this->xmx; }

//! get ymx attribute value
/*!
\return The ymx attribute value - float type
*/
float Hydyn::getYmx() { return this->ymx; }

//! get zmx attribute value
/*!
\return The zmx attribute value - float type
*/
float Hydyn::getZmx() { return this->zmx; }

//! get mmx attribute value
/*!
\return The mmx attribute value - int type
*/
int Hydyn::getMmx() { return this->mmx; }

//! get mmy attribute value
/*!
\return The mmy attribute value - int type
*/
int Hydyn::getMmy() { return this->mmy; }

//! get mmz attribute value
/*!
\return The mmz attribute value - int type
*/
int Hydyn::getMmz() { return this->mmz; }

//! get ratioeex attribute value
/*!
\return The ratioeex attribute value - Array<float> * type
*/
Array<float> *Hydyn::getRatioeex() { return this->ratioeex; }

//! get velc attribute value
/*!
\return The velc attribute value - Array<double> * type
*/
Array<double> *Hydyn::getVelc() { return this->velc; }

//! get epsc attribute value
/*!
\return The epsc attribute value - Array<double> * type
*/
Array<double> *Hydyn::getEpsc() { return this->epsc; }

//! get sigc attribute value
/*!
\return The sigc attribute value - Array<double> * type
*/
Array<double> *Hydyn::getSigc() { return this->sigc; }

//! get barc attribute value
/*!
\return The barc attribute value - Array<double> * type
*/
Array<double> *Hydyn::getBarc() { return this->barc; }

// c++ attribute mutator methods
//! set nxhy attribute value
/*!
\param nxhyParam - the nxhy attribute value - int type
*/
void Hydyn::setNxhy(int nxhyParam) {
  this->nxhy = nxhyParam;
  // set value from fortran common block
  this->setFortranNxhy();
}

//! set nyhy attribute value
/*!
\param nyhyParam - the nyhy attribute value - int type
*/
void Hydyn::setNyhy(int nyhyParam) {
  this->nyhy = nyhyParam;
  // set value from fortran common block
  this->setFortranNyhy();
}

//! set nzhy attribute value
/*!
\param nzhyParam - the nzhy attribute value - int type
*/
void Hydyn::setNzhy(int nzhyParam) {
  this->nzhy = nzhyParam;
  // set value from fortran common block
  this->setFortranNzhy();
}

//! set ntauhy attribute value
/*!
\param ntauhyParam - the ntauhy attribute value - int type
*/
void Hydyn::setNtauhy(int ntauhyParam) {
  this->ntauhy = ntauhyParam;
  // set value from fortran common block
  this->setFortranNtauhy();
}

//! set dtauhy attribute value
/*!
\param dtauhyParam - the dtauhy attribute value - float type
*/
void Hydyn::setDtauhy(float dtauhyParam) {
  this->dtauhy = dtauhyParam;
  // set value from fortran common block
  this->setFortranDtauhy();
}

//! set xminhy attribute value
/*!
\param xminhyParam - the xminhy attribute value - float type
*/
void Hydyn::setXminhy(float xminhyParam) {
  this->xminhy = xminhyParam;
  // set value from fortran common block
  this->setFortranXminhy();
}

//! set xmaxhy attribute value
/*!
\param xmaxhyParam - the xmaxhy attribute value - float type
*/
void Hydyn::setXmaxhy(float xmaxhyParam) {
  this->xmaxhy = xmaxhyParam;
  // set value from fortran common block
  this->setFortranXmaxhy();
}

//! set yminhy attribute value
/*!
\param yminhyParam - the yminhy attribute value - float type
*/
void Hydyn::setYminhy(float yminhyParam) {
  this->yminhy = yminhyParam;
  // set value from fortran common block
  this->setFortranYminhy();
}

//! set ymaxhy attribute value
/*!
\param ymaxhyParam - the ymaxhy attribute value - float type
*/
void Hydyn::setYmaxhy(float ymaxhyParam) {
  this->ymaxhy = ymaxhyParam;
  // set value from fortran common block
  this->setFortranYmaxhy();
}

//! set zminhy attribute value
/*!
\param zminhyParam - the zminhy attribute value - float type
*/
void Hydyn::setZminhy(float zminhyParam) {
  this->zminhy = zminhyParam;
  // set value from fortran common block
  this->setFortranZminhy();
}

//! set zmaxhy attribute value
/*!
\param zmaxhyParam - the zmaxhy attribute value - float type
*/
void Hydyn::setZmaxhy(float zmaxhyParam) {
  this->zmaxhy = zmaxhyParam;
  // set value from fortran common block
  this->setFortranZmaxhy();
}

//! set tfrout attribute value
/*!
\param tfroutParam - the tfrout attribute value - float type
*/
void Hydyn::setTfrout(float tfroutParam) {
  this->tfrout = tfroutParam;
  // set value from fortran common block
  this->setFortranTfrout();
}

//! set ioeos attribute value
/*!
\param ioeosParam - the ioeos attribute value - int type
*/
void Hydyn::setIoeos(int ioeosParam) {
  this->ioeos = ioeosParam;
  // set value from fortran common block
  this->setFortranIoeos();
}

//! set iozerof attribute value
/*!
\param iozerofParam - the iozerof attribute value - int type
*/
void Hydyn::setIozerof(int iozerofParam) {
  this->iozerof = iozerofParam;
  // set value from fortran common block
  this->setFortranIozerof();
}

//! set ihlle attribute value
/*!
\param ihlleParam - the ihlle attribute value - int type
*/
void Hydyn::setIhlle(int ihlleParam) {
  this->ihlle = ihlleParam;
  // set value from fortran common block
  this->setFortranIhlle();
}

//! set jhlle attribute value
/*!
\param jhlleParam - the jhlle attribute value - int type
*/
void Hydyn::setJhlle(int jhlleParam) {
  this->jhlle = jhlleParam;
  // set value from fortran common block
  this->setFortranJhlle();
}

//! set ntaumx attribute value
/*!
\param ntaumxParam - the ntaumx attribute value - int type
*/
void Hydyn::setNtaumx(int ntaumxParam) {
  this->ntaumx = ntaumxParam;
  // set value from fortran common block
  this->setFortranNtaumx();
}

//! set tfo attribute value
/*!
\param tfoParam - the tfo attribute value - float type
*/
void Hydyn::setTfo(float tfoParam) {
  this->tfo = tfoParam;
  // set value from fortran common block
  this->setFortranTfo();
}

//! set ienvar attribute value
/*!
\param ienvarParam - the ienvar attribute value - int type
*/
void Hydyn::setIenvar(int ienvarParam) {
  this->ienvar = ienvarParam;
  // set value from fortran common block
  this->setFortranIenvar();
}

//! set ifaahlle attribute value
/*!
\param ifaahlleParam - the ifaahlle attribute value - int type
*/
void Hydyn::setIfaahlle(int ifaahlleParam) {
  this->ifaahlle = ifaahlleParam;
  // set value from fortran common block
  this->setFortranIfaahlle();
}

//! set ifathlle attribute value
/*!
\param ifathlleParam - the ifathlle attribute value - int type
*/
void Hydyn::setIfathlle(int ifathlleParam) {
  this->ifathlle = ifathlleParam;
  // set value from fortran common block
  this->setFortranIfathlle();
}

//! set ifazhlle attribute value
/*!
\param ifazhlleParam - the ifazhlle attribute value - int type
*/
void Hydyn::setIfazhlle(int ifazhlleParam) {
  this->ifazhlle = ifazhlleParam;
  // set value from fortran common block
  this->setFortranIfazhlle();
}

//! set epsfin attribute value
/*!
\param epsfinParam - the epsfin attribute value - float type
*/
void Hydyn::setEpsfin(float epsfinParam) {
  this->epsfin = epsfinParam;
  // set value from fortran common block
  this->setFortranEpsfin();
}

//! set tauhy attribute value
/*!
\param tauhyParam - the tauhy attribute value - Array<float> * type
*/
void Hydyn::setTauhy(Array<float> *tauhyParam) { this->tauhy = tauhyParam; }

//! set eetau attribute value
/*!
\param eetauParam - the eetau attribute value - float type
*/
void Hydyn::setEetau(float eetauParam) { this->eetau = eetauParam; }

//! set eetau2 attribute value
/*!
\param eetau2Param - the eetau2 attribute value - float type
*/
void Hydyn::setEetau2(float eetau2Param) { this->eetau2 = eetau2Param; }

//! set emx attribute value
/*!
\param emxParam - the emx attribute value - float type
*/
void Hydyn::setEmx(float emxParam) { this->emx = emxParam; }

//! set xmx attribute value
/*!
\param xmxParam - the xmx attribute value - float type
*/
void Hydyn::setXmx(float xmxParam) { this->xmx = xmxParam; }

//! set ymx attribute value
/*!
\param ymxParam - the ymx attribute value - float type
*/
void Hydyn::setYmx(float ymxParam) { this->ymx = ymxParam; }

//! set zmx attribute value
/*!
\param zmxParam - the zmx attribute value - float type
*/
void Hydyn::setZmx(float zmxParam) { this->zmx = zmxParam; }

//! set mmx attribute value
/*!
\param mmxParam - the mmx attribute value - int type
*/
void Hydyn::setMmx(int mmxParam) { this->mmx = mmxParam; }

//! set mmy attribute value
/*!
\param mmyParam - the mmy attribute value - int type
*/
void Hydyn::setMmy(int mmyParam) { this->mmy = mmyParam; }

//! set mmz attribute value
/*!
\param mmzParam - the mmz attribute value - int type
*/
void Hydyn::setMmz(int mmzParam) { this->mmz = mmzParam; }

//! set ratioeex attribute value
/*!
\param ratioeexParam - the ratioeex attribute value - Array<float> * type
*/
void Hydyn::setRatioeex(Array<float> *ratioeexParam) { this->ratioeex = ratioeexParam; }

//! set velc attribute value
/*!
\param velcParam - the velc attribute value - Array<double> * type
*/
void Hydyn::setVelc(Array<double> *velcParam) { this->velc = velcParam; }

//! set epsc attribute value
/*!
\param epscParam - the epsc attribute value - Array<double> * type
*/
void Hydyn::setEpsc(Array<double> *epscParam) { this->epsc = epscParam; }

//! set sigc attribute value
/*!
\param sigcParam - the sigc attribute value - Array<double> * type
*/
void Hydyn::setSigc(Array<double> *sigcParam) { this->sigc = sigcParam; }

//! set barc attribute value
/*!
\param barcParam - the barc attribute value - Array<double> * type
*/
void Hydyn::setBarc(Array<double> *barcParam) { this->barc = barcParam; }

// c++ array attribute methods
//! create tauhy attribute array of float type
/*!
\param ntau - indice - int type
*/
void Hydyn::createTauhy(int ntau) {
  int dim = 1;
  int shape[1] = {ntau};
  this->tauhy = new Array<float>(dim, shape);
}

//! get tauhy attribute array value
/*!
\param ntau - indice - int type
\return The tauhy attribute value - float type
*/
float Hydyn::getTauhyValue(int ntau){
  int index[1] = {ntau};
  return this->tauhy->get(index);
}

//! set tauhy attribute array value
/*!
\param ntau - indice - int type
\param value - the value to be set - float type
*/
void Hydyn::setTauhyValue(int ntau, float value) {
  int index[1] = {ntau};
  this->tauhy->set(index, value);
}

//! create tauhy initialize method
void Hydyn::initializeTauhy(float value) {
  this->tauhy->initialize(value);
}

//! create ratioeex attribute array of float type
/*!
\param ntau - indice - int type
*/
void Hydyn::createRatioeex(int ntau) {
  int dim = 1;
  int shape[1] = {ntau};
  this->ratioeex = new Array<float>(dim, shape);
}

//! get ratioeex attribute array value
/*!
\param ntau - indice - int type
\return The ratioeex attribute value - float type
*/
float Hydyn::getRatioeexValue(int ntau){
  int index[1] = {ntau};
  return this->ratioeex->get(index);
}

//! set ratioeex attribute array value
/*!
\param ntau - indice - int type
\param value - the value to be set - float type
*/
void Hydyn::setRatioeexValue(int ntau, float value) {
  int index[1] = {ntau};
  this->ratioeex->set(index, value);
}

//! create ratioeex initialize method
void Hydyn::initializeRatioeex(float value) {
  this->ratioeex->initialize(value);
}

//! create velc attribute array of double type
/*!
\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/
void Hydyn::createVelc(int n, int neta, int ntau, int nx, int ny) {
  int dim = 5;
  int shape[5] = {n, neta, ntau, nx, ny};
  this->velc = new Array<double>(dim, shape);
}

//! get velc attribute array value
/*!
\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\return The velc attribute value - double type
*/
double Hydyn::getVelcValue(int n, int neta, int ntau, int nx, int ny){
  int index[5] = {n, neta, ntau, nx, ny};
  return this->velc->get(index);
}

//! set velc attribute array value
/*!
\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\param value - the value to be set - double type
*/
void Hydyn::setVelcValue(int n, int neta, int ntau, int nx, int ny, double value) {
  int index[5] = {n, neta, ntau, nx, ny};
  this->velc->set(index, value);
}

//! create velc initialize method
void Hydyn::initializeVelc(double value) {
  this->velc->initialize(value);
}

//! create epsc attribute array of double type
/*!
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/
void Hydyn::createEpsc(int neta, int ntau, int nx, int ny) {
  int dim = 4;
  int shape[4] = {neta, ntau, nx, ny};
  this->epsc = new Array<double>(dim, shape);
}

//! get epsc attribute array value
/*!
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\return The epsc attribute value - double type
*/
double Hydyn::getEpscValue(int neta, int ntau, int nx, int ny){
  int index[4] = {neta, ntau, nx, ny};
  return this->epsc->get(index);
}

//! set epsc attribute array value
/*!
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\param value - the value to be set - double type
*/
void Hydyn::setEpscValue(int neta, int ntau, int nx, int ny, double value) {
  int index[4] = {neta, ntau, nx, ny};
  this->epsc->set(index, value);
}

//! create epsc initialize method
void Hydyn::initializeEpsc(double value) {
  this->epsc->initialize(value);
}

//! create sigc attribute array of double type
/*!
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/
void Hydyn::createSigc(int neta, int ntau, int nx, int ny) {
  int dim = 4;
  int shape[4] = {neta, ntau, nx, ny};
  this->sigc = new Array<double>(dim, shape);
}

//! get sigc attribute array value
/*!
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\return The sigc attribute value - double type
*/
double Hydyn::getSigcValue(int neta, int ntau, int nx, int ny){
  int index[4] = {neta, ntau, nx, ny};
  return this->sigc->get(index);
}

//! set sigc attribute array value
/*!
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\param value - the value to be set - double type
*/
void Hydyn::setSigcValue(int neta, int ntau, int nx, int ny, double value) {
  int index[4] = {neta, ntau, nx, ny};
  this->sigc->set(index, value);
}

//! create sigc initialize method
void Hydyn::initializeSigc(double value) {
  this->sigc->initialize(value);
}

//! create barc attribute array of double type
/*!
\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/
void Hydyn::createBarc(int n, int neta, int ntau, int nx, int ny) {
  int dim = 5;
  int shape[5] = {n, neta, ntau, nx, ny};
  this->barc = new Array<double>(dim, shape);
}

//! get barc attribute array value
/*!
\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\return The barc attribute value - double type
*/
double Hydyn::getBarcValue(int n, int neta, int ntau, int nx, int ny){
  int index[5] = {n, neta, ntau, nx, ny};
  return this->barc->get(index);
}

//! set barc attribute array value
/*!
\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
\param value - the value to be set - double type
*/
void Hydyn::setBarcValue(int n, int neta, int ntau, int nx, int ny, double value) {
  int index[5] = {n, neta, ntau, nx, ny};
  this->barc->set(index, value);
}

//! create barc initialize method
void Hydyn::initializeBarc(double value) {
  this->barc->initialize(value);
}

//! get nxhy, nyhy, nzhy, ntauhy attribute values
/*!
\param nxhyParam - the nxhy attribute address - int pointer type
\param nyhyParam - the nyhy attribute address - int pointer type
\param nzhyParam - the nzhy attribute address - int pointer type
\param ntauhyParam - the ntauhy attribute address - int pointer type
*/
void Hydyn::getDim(int *nxhyParam, int *nyhyParam, int *nzhyParam, int *ntauhyParam) {
  *nxhyParam = this->getNxhy();
  *nyhyParam = this->getNyhy();
  *nzhyParam = this->getNzhy();
  *ntauhyParam = this->getNtauhy();
}

//! set nxhy, nyhy, nzhy, ntauhy attribute values
/*!
\param nxhyParam - the nxhy attribute address - int type
\param nyhyParam - the nyhy attribute address - int type
\param nzhyParam - the nzhy attribute address - int type
\param ntauhyParam - the ntauhy attribute address - int type
*/
void Hydyn::setDim(int nxhyParam, int nyhyParam, int nzhyParam, int ntauhyParam) {
  this->setNxhy(nxhyParam);
  this->setNyhy(nyhyParam);
  this->setNzhy(nzhyParam);
  this->setNtauhy(ntauhyParam);
}

//! get xminhy, xmaxhy, yminhy, ymaxhy, zminhy, zmaxhy attribute values
/*!
\param xminhyParam - the xminhy attribute address - float pointer type
\param xmaxhyParam - the xmaxhy attribute address - float pointer type
\param yminhyParam - the yminhy attribute address - float pointer type
\param ymaxhyParam - the ymaxhy attribute address - float pointer type
\param zminhyParam - the zminhy attribute address - float pointer type
\param zmaxhyParam - the zmaxhy attribute address - float pointer type
*/
void Hydyn::getBound(float *xminhyParam, float *xmaxhyParam, float *yminhyParam, float *ymaxhyParam, float *zminhyParam, float *zmaxhyParam) {
  *xminhyParam = this->getXminhy();
  *xmaxhyParam = this->getXmaxhy();
  *yminhyParam = this->getYminhy();
  *ymaxhyParam = this->getYmaxhy();
  *zminhyParam = this->getZminhy();
  *zmaxhyParam = this->getZmaxhy();
}

//! set xminhy, xmaxhy, yminhy, ymaxhy, zminhy, zmaxhy attribute values
/*!
\param xminhyParam - the xminhy attribute address - float type
\param xmaxhyParam - the xmaxhy attribute address - float type
\param yminhyParam - the yminhy attribute address - float type
\param ymaxhyParam - the ymaxhy attribute address - float type
\param zminhyParam - the zminhy attribute address - float type
\param zmaxhyParam - the zmaxhy attribute address - float type
*/
void Hydyn::setBound(float xminhyParam, float xmaxhyParam, float yminhyParam, float ymaxhyParam, float zminhyParam, float zmaxhyParam) {
  this->setXminhy(xminhyParam);
  this->setXmaxhy(xmaxhyParam);
  this->setYminhy(yminhyParam);
  this->setYmaxhy(ymaxhyParam);
  this->setZminhy(zminhyParam);
  this->setZmaxhy(zmaxhyParam);
}

//! get ioeos, iozerof attribute values
/*!
\param ioeosParam - the ioeos attribute address - int pointer type
\param iozerofParam - the iozerof attribute address - int pointer type
*/
void Hydyn::getIo(int *ioeosParam, int *iozerofParam) {
  *ioeosParam = this->getIoeos();
  *iozerofParam = this->getIozerof();
}

//! set ioeos, iozerof attribute values
/*!
\param ioeosParam - the ioeos attribute address - int type
\param iozerofParam - the iozerof attribute address - int type
*/
void Hydyn::setIo(int ioeosParam, int iozerofParam) {
  this->setIoeos(ioeosParam);
  this->setIozerof(iozerofParam);
}

// Hydyn instance pointer
Hydyn *hydyn;

// Fortran / C wrapper functions
extern "C" {
//! create Hydyn instance
void createhydyn_() { hydyn = new Hydyn(); }

//! delete Hydyn instance
void destroyhydyn_() { delete hydyn; }

//! get nxhy attribute value
/*!
\return The nxhy attribute value - int type
*/
int gethydynnxhy_() { return hydyn->getNxhy(); }

//! set nxhy attribute value
/*!
\param nxhyParam - the nxhy attribute value - int type
*/
void sethydynnxhy_(int *nxhyParam) { hydyn->setNxhy(*nxhyParam); }

//! get nyhy attribute value
/*!
\return The nyhy attribute value - int type
*/
int gethydynnyhy_() { return hydyn->getNyhy(); }

//! set nyhy attribute value
/*!
\param nyhyParam - the nyhy attribute value - int type
*/
void sethydynnyhy_(int *nyhyParam) { hydyn->setNyhy(*nyhyParam); }

//! get nzhy attribute value
/*!
\return The nzhy attribute value - int type
*/
int gethydynnzhy_() { return hydyn->getNzhy(); }

//! set nzhy attribute value
/*!
\param nzhyParam - the nzhy attribute value - int type
*/
void sethydynnzhy_(int *nzhyParam) { hydyn->setNzhy(*nzhyParam); }

//! get ntauhy attribute value
/*!
\return The ntauhy attribute value - int type
*/
int gethydynntauhy_() { return hydyn->getNtauhy(); }

//! set ntauhy attribute value
/*!
\param ntauhyParam - the ntauhy attribute value - int type
*/
void sethydynntauhy_(int *ntauhyParam) { hydyn->setNtauhy(*ntauhyParam); }

//! get dtauhy attribute value
/*!
\return The dtauhy attribute value - float type
*/
float gethydyndtauhy_() { return hydyn->getDtauhy(); }

//! set dtauhy attribute value
/*!
\param dtauhyParam - the dtauhy attribute value - float type
*/
void sethydyndtauhy_(float *dtauhyParam) { hydyn->setDtauhy(*dtauhyParam); }

//! get xminhy attribute value
/*!
\return The xminhy attribute value - float type
*/
float gethydynxminhy_() { return hydyn->getXminhy(); }

//! set xminhy attribute value
/*!
\param xminhyParam - the xminhy attribute value - float type
*/
void sethydynxminhy_(float *xminhyParam) { hydyn->setXminhy(*xminhyParam); }

//! get xmaxhy attribute value
/*!
\return The xmaxhy attribute value - float type
*/
float gethydynxmaxhy_() { return hydyn->getXmaxhy(); }

//! set xmaxhy attribute value
/*!
\param xmaxhyParam - the xmaxhy attribute value - float type
*/
void sethydynxmaxhy_(float *xmaxhyParam) { hydyn->setXmaxhy(*xmaxhyParam); }

//! get yminhy attribute value
/*!
\return The yminhy attribute value - float type
*/
float gethydynyminhy_() { return hydyn->getYminhy(); }

//! set yminhy attribute value
/*!
\param yminhyParam - the yminhy attribute value - float type
*/
void sethydynyminhy_(float *yminhyParam) { hydyn->setYminhy(*yminhyParam); }

//! get ymaxhy attribute value
/*!
\return The ymaxhy attribute value - float type
*/
float gethydynymaxhy_() { return hydyn->getYmaxhy(); }

//! set ymaxhy attribute value
/*!
\param ymaxhyParam - the ymaxhy attribute value - float type
*/
void sethydynymaxhy_(float *ymaxhyParam) { hydyn->setYmaxhy(*ymaxhyParam); }

//! get zminhy attribute value
/*!
\return The zminhy attribute value - float type
*/
float gethydynzminhy_() { return hydyn->getZminhy(); }

//! set zminhy attribute value
/*!
\param zminhyParam - the zminhy attribute value - float type
*/
void sethydynzminhy_(float *zminhyParam) { hydyn->setZminhy(*zminhyParam); }

//! get zmaxhy attribute value
/*!
\return The zmaxhy attribute value - float type
*/
float gethydynzmaxhy_() { return hydyn->getZmaxhy(); }

//! set zmaxhy attribute value
/*!
\param zmaxhyParam - the zmaxhy attribute value - float type
*/
void sethydynzmaxhy_(float *zmaxhyParam) { hydyn->setZmaxhy(*zmaxhyParam); }

//! get tfrout attribute value
/*!
\return The tfrout attribute value - float type
*/
float gethydyntfrout_() { return hydyn->getTfrout(); }

//! set tfrout attribute value
/*!
\param tfroutParam - the tfrout attribute value - float type
*/
void sethydyntfrout_(float *tfroutParam) { hydyn->setTfrout(*tfroutParam); }

//! get ioeos attribute value
/*!
\return The ioeos attribute value - int type
*/
int gethydynioeos_() { return hydyn->getIoeos(); }

//! set ioeos attribute value
/*!
\param ioeosParam - the ioeos attribute value - int type
*/
void sethydynioeos_(int *ioeosParam) { hydyn->setIoeos(*ioeosParam); }

//! get iozerof attribute value
/*!
\return The iozerof attribute value - int type
*/
int gethydyniozerof_() { return hydyn->getIozerof(); }

//! set iozerof attribute value
/*!
\param iozerofParam - the iozerof attribute value - int type
*/
void sethydyniozerof_(int *iozerofParam) { hydyn->setIozerof(*iozerofParam); }

//! get ihlle attribute value
/*!
\return The ihlle attribute value - int type
*/
int gethydynihlle_() { return hydyn->getIhlle(); }

//! set ihlle attribute value
/*!
\param ihlleParam - the ihlle attribute value - int type
*/
void sethydynihlle_(int *ihlleParam) { hydyn->setIhlle(*ihlleParam); }

//! get jhlle attribute value
/*!
\return The jhlle attribute value - int type
*/
int gethydynjhlle_() { return hydyn->getJhlle(); }

//! set jhlle attribute value
/*!
\param jhlleParam - the jhlle attribute value - int type
*/
void sethydynjhlle_(int *jhlleParam) { hydyn->setJhlle(*jhlleParam); }

//! get ntaumx attribute value
/*!
\return The ntaumx attribute value - int type
*/
int gethydynntaumx_() { return hydyn->getNtaumx(); }

//! set ntaumx attribute value
/*!
\param ntaumxParam - the ntaumx attribute value - int type
*/
void sethydynntaumx_(int *ntaumxParam) { hydyn->setNtaumx(*ntaumxParam); }

//! get tfo attribute value
/*!
\return The tfo attribute value - float type
*/
float gethydyntfo_() { return hydyn->getTfo(); }

//! set tfo attribute value
/*!
\param tfoParam - the tfo attribute value - float type
*/
void sethydyntfo_(float *tfoParam) { hydyn->setTfo(*tfoParam); }

//! get ienvar attribute value
/*!
\return The ienvar attribute value - int type
*/
int gethydynienvar_() { return hydyn->getIenvar(); }

//! set ienvar attribute value
/*!
\param ienvarParam - the ienvar attribute value - int type
*/
void sethydynienvar_(int *ienvarParam) { hydyn->setIenvar(*ienvarParam); }

//! get ifaahlle attribute value
/*!
\return The ifaahlle attribute value - int type
*/
int gethydynifaahlle_() { return hydyn->getIfaahlle(); }

//! set ifaahlle attribute value
/*!
\param ifaahlleParam - the ifaahlle attribute value - int type
*/
void sethydynifaahlle_(int *ifaahlleParam) { hydyn->setIfaahlle(*ifaahlleParam); }

//! get ifathlle attribute value
/*!
\return The ifathlle attribute value - int type
*/
int gethydynifathlle_() { return hydyn->getIfathlle(); }

//! set ifathlle attribute value
/*!
\param ifathlleParam - the ifathlle attribute value - int type
*/
void sethydynifathlle_(int *ifathlleParam) { hydyn->setIfathlle(*ifathlleParam); }

//! get ifazhlle attribute value
/*!
\return The ifazhlle attribute value - int type
*/
int gethydynifazhlle_() { return hydyn->getIfazhlle(); }

//! set ifazhlle attribute value
/*!
\param ifazhlleParam - the ifazhlle attribute value - int type
*/
void sethydynifazhlle_(int *ifazhlleParam) { hydyn->setIfazhlle(*ifazhlleParam); }

//! get epsfin attribute value
/*!
\return The epsfin attribute value - float type
*/
float gethydynepsfin_() { return hydyn->getEpsfin(); }

//! set epsfin attribute value
/*!
\param epsfinParam - the epsfin attribute value - float type
*/
void sethydynepsfin_(float *epsfinParam) { hydyn->setEpsfin(*epsfinParam); }

//! create tauhy array attribute
/*!\param ntau - indice - int type
*/void createhydyntauhy_(int* ntau) {
  hydyn->createTauhy(*ntau);
}

//! initialize tauhy array attribute
void initializehydyntauhy_(float *value) {
  hydyn->initializeTauhy(*value);
}

//! set tauhy array value
/*!
\param ntau - indice - int pointer type
\param value - value to set - float pointer type
*/
void sethydyntauhy_(int *ntau, float *value) {
  hydyn->setTauhyValue(*ntau, *value);
}

//! get tauhy array value
/*!
\param ntau indice - int pointer type
\return value - value to get - float pointer type
*/
float gethydyntauhy_(int *ntau) {
  return hydyn->getTauhyValue(*ntau);
}

//! delete tauhy array
void destroyhydyntauhy_() { delete hydyn->getTauhy(); }

//! get eetau attribute value
/*!
\return The eetau attribute value - float type
*/
float gethydyneetau_() { return hydyn->getEetau(); }

//! set eetau attribute value
/*!
\param eetauParam - the eetau attribute value - float type
*/
void sethydyneetau_(float *eetauParam) { hydyn->setEetau(*eetauParam); }

//! get eetau2 attribute value
/*!
\return The eetau2 attribute value - float type
*/
float gethydyneetau2_() { return hydyn->getEetau2(); }

//! set eetau2 attribute value
/*!
\param eetau2Param - the eetau2 attribute value - float type
*/
void sethydyneetau2_(float *eetau2Param) { hydyn->setEetau2(*eetau2Param); }

//! get emx attribute value
/*!
\return The emx attribute value - float type
*/
float gethydynemx_() { return hydyn->getEmx(); }

//! set emx attribute value
/*!
\param emxParam - the emx attribute value - float type
*/
void sethydynemx_(float *emxParam) { hydyn->setEmx(*emxParam); }

//! get xmx attribute value
/*!
\return The xmx attribute value - float type
*/
float gethydynxmx_() { return hydyn->getXmx(); }

//! set xmx attribute value
/*!
\param xmxParam - the xmx attribute value - float type
*/
void sethydynxmx_(float *xmxParam) { hydyn->setXmx(*xmxParam); }

//! get ymx attribute value
/*!
\return The ymx attribute value - float type
*/
float gethydynymx_() { return hydyn->getYmx(); }

//! set ymx attribute value
/*!
\param ymxParam - the ymx attribute value - float type
*/
void sethydynymx_(float *ymxParam) { hydyn->setYmx(*ymxParam); }

//! get zmx attribute value
/*!
\return The zmx attribute value - float type
*/
float gethydynzmx_() { return hydyn->getZmx(); }

//! set zmx attribute value
/*!
\param zmxParam - the zmx attribute value - float type
*/
void sethydynzmx_(float *zmxParam) { hydyn->setZmx(*zmxParam); }

//! get mmx attribute value
/*!
\return The mmx attribute value - int type
*/
int gethydynmmx_() { return hydyn->getMmx(); }

//! set mmx attribute value
/*!
\param mmxParam - the mmx attribute value - int type
*/
void sethydynmmx_(int *mmxParam) { hydyn->setMmx(*mmxParam); }

//! get mmy attribute value
/*!
\return The mmy attribute value - int type
*/
int gethydynmmy_() { return hydyn->getMmy(); }

//! set mmy attribute value
/*!
\param mmyParam - the mmy attribute value - int type
*/
void sethydynmmy_(int *mmyParam) { hydyn->setMmy(*mmyParam); }

//! get mmz attribute value
/*!
\return The mmz attribute value - int type
*/
int gethydynmmz_() { return hydyn->getMmz(); }

//! set mmz attribute value
/*!
\param mmzParam - the mmz attribute value - int type
*/
void sethydynmmz_(int *mmzParam) { hydyn->setMmz(*mmzParam); }

//! create ratioeex array attribute
/*!\param ntau - indice - int type
*/void createhydynratioeex_(int* ntau) {
  hydyn->createRatioeex(*ntau);
}

//! initialize ratioeex array attribute
void initializehydynratioeex_(float *value) {
  hydyn->initializeRatioeex(*value);
}

//! set ratioeex array value
/*!
\param ntau - indice - int pointer type
\param value - value to set - float pointer type
*/
void sethydynratioeex_(int *ntau, float *value) {
  hydyn->setRatioeexValue(*ntau, *value);
}

//! get ratioeex array value
/*!
\param ntau indice - int pointer type
\return value - value to get - float pointer type
*/
float gethydynratioeex_(int *ntau) {
  return hydyn->getRatioeexValue(*ntau);
}

//! delete ratioeex array
void destroyhydynratioeex_() { delete hydyn->getRatioeex(); }

//! create velc array attribute
/*!\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/void createhydynvelc_(int* n, int* neta, int* ntau, int* nx, int* ny) {
  hydyn->createVelc(*n, *neta, *ntau, *nx, *ny);
}

//! initialize velc array attribute
void initializehydynvelc_(double *value) {
  hydyn->initializeVelc(*value);
}

//! set velc array value
/*!
\param n - indice - int pointer type
\param neta - indice - int pointer type
\param ntau - indice - int pointer type
\param nx - indice - int pointer type
\param ny - indice - int pointer type
\param value - value to set - double pointer type
*/
void sethydynvelc_(int *n, int *neta, int *ntau, int *nx, int *ny, double *value) {
  hydyn->setVelcValue(*n, *neta, *ntau, *nx, *ny, *value);
}

//! get velc array value
/*!
\param n indice - int pointer type
\param neta indice - int pointer type
\param ntau indice - int pointer type
\param nx indice - int pointer type
\param ny indice - int pointer type
\return value - value to get - double pointer type
*/
double gethydynvelc_(int *n, int *neta, int *ntau, int *nx, int *ny) {
  return hydyn->getVelcValue(*n, *neta, *ntau, *nx, *ny);
}

//! delete velc array
void destroyhydynvelc_() { delete hydyn->getVelc(); }

//! create epsc array attribute
/*!\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/void createhydynepsc_(int* neta, int* ntau, int* nx, int* ny) {
  hydyn->createEpsc(*neta, *ntau, *nx, *ny);
}

//! initialize epsc array attribute
void initializehydynepsc_(double *value) {
  hydyn->initializeEpsc(*value);
}

//! set epsc array value
/*!
\param neta - indice - int pointer type
\param ntau - indice - int pointer type
\param nx - indice - int pointer type
\param ny - indice - int pointer type
\param value - value to set - double pointer type
*/
void sethydynepsc_(int *neta, int *ntau, int *nx, int *ny, double *value) {
  hydyn->setEpscValue(*neta, *ntau, *nx, *ny, *value);
}

//! get epsc array value
/*!
\param neta indice - int pointer type
\param ntau indice - int pointer type
\param nx indice - int pointer type
\param ny indice - int pointer type
\return value - value to get - double pointer type
*/
double gethydynepsc_(int *neta, int *ntau, int *nx, int *ny) {
  return hydyn->getEpscValue(*neta, *ntau, *nx, *ny);
}

//! delete epsc array
void destroyhydynepsc_() { delete hydyn->getEpsc(); }

//! create sigc array attribute
/*!\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/void createhydynsigc_(int* neta, int* ntau, int* nx, int* ny) {
  hydyn->createSigc(*neta, *ntau, *nx, *ny);
}

//! initialize sigc array attribute
void initializehydynsigc_(double *value) {
  hydyn->initializeSigc(*value);
}

//! set sigc array value
/*!
\param neta - indice - int pointer type
\param ntau - indice - int pointer type
\param nx - indice - int pointer type
\param ny - indice - int pointer type
\param value - value to set - double pointer type
*/
void sethydynsigc_(int *neta, int *ntau, int *nx, int *ny, double *value) {
  hydyn->setSigcValue(*neta, *ntau, *nx, *ny, *value);
}

//! get sigc array value
/*!
\param neta indice - int pointer type
\param ntau indice - int pointer type
\param nx indice - int pointer type
\param ny indice - int pointer type
\return value - value to get - double pointer type
*/
double gethydynsigc_(int *neta, int *ntau, int *nx, int *ny) {
  return hydyn->getSigcValue(*neta, *ntau, *nx, *ny);
}

//! delete sigc array
void destroyhydynsigc_() { delete hydyn->getSigc(); }

//! create barc array attribute
/*!\param n - indice - int type
\param neta - indice - int type
\param ntau - indice - int type
\param nx - indice - int type
\param ny - indice - int type
*/void createhydynbarc_(int* n, int* neta, int* ntau, int* nx, int* ny) {
  hydyn->createBarc(*n, *neta, *ntau, *nx, *ny);
}

//! initialize barc array attribute
void initializehydynbarc_(double *value) {
  hydyn->initializeBarc(*value);
}

//! set barc array value
/*!
\param n - indice - int pointer type
\param neta - indice - int pointer type
\param ntau - indice - int pointer type
\param nx - indice - int pointer type
\param ny - indice - int pointer type
\param value - value to set - double pointer type
*/
void sethydynbarc_(int *n, int *neta, int *ntau, int *nx, int *ny, double *value) {
  hydyn->setBarcValue(*n, *neta, *ntau, *nx, *ny, *value);
}

//! get barc array value
/*!
\param n indice - int pointer type
\param neta indice - int pointer type
\param ntau indice - int pointer type
\param nx indice - int pointer type
\param ny indice - int pointer type
\return value - value to get - double pointer type
*/
double gethydynbarc_(int *n, int *neta, int *ntau, int *nx, int *ny) {
  return hydyn->getBarcValue(*n, *neta, *ntau, *nx, *ny);
}

//! delete barc array
void destroyhydynbarc_() { delete hydyn->getBarc(); }

//! get nxhy, nyhy, nzhy, ntauhy attribute values
/*!
\param nxhyParam - the nxhy attribute address - int pointer type
\param nyhyParam - the nyhy attribute address - int pointer type
\param nzhyParam - the nzhy attribute address - int pointer type
\param ntauhyParam - the ntauhy attribute address - int pointer type
*/
void gethydyndim_(int* nxhyParam, int* nyhyParam, int* nzhyParam, int* ntauhyParam) {
  hydyn->getDim(nxhyParam,nyhyParam,nzhyParam,ntauhyParam);
}

//! set nxhy, nyhy, nzhy, ntauhy attribute values
/*!
\param nxhyParam - the nxhy attribute address - int pointer type
\param nyhyParam - the nyhy attribute address - int pointer type
\param nzhyParam - the nzhy attribute address - int pointer type
\param ntauhyParam - the ntauhy attribute address - int pointer type
*/
void sethydyndim_(int* nxhyParam, int* nyhyParam, int* nzhyParam, int* ntauhyParam) {
  hydyn->setDim(*nxhyParam,*nyhyParam,*nzhyParam,*ntauhyParam);
}

//! get xminhy, xmaxhy, yminhy, ymaxhy, zminhy, zmaxhy attribute values
/*!
\param xminhyParam - the xminhy attribute address - float pointer type
\param xmaxhyParam - the xmaxhy attribute address - float pointer type
\param yminhyParam - the yminhy attribute address - float pointer type
\param ymaxhyParam - the ymaxhy attribute address - float pointer type
\param zminhyParam - the zminhy attribute address - float pointer type
\param zmaxhyParam - the zmaxhy attribute address - float pointer type
*/
void gethydynbound_(float* xminhyParam, float* xmaxhyParam, float* yminhyParam, float* ymaxhyParam, float* zminhyParam, float* zmaxhyParam) {
  hydyn->getBound(xminhyParam,xmaxhyParam,yminhyParam,ymaxhyParam,zminhyParam,zmaxhyParam);
}

//! set xminhy, xmaxhy, yminhy, ymaxhy, zminhy, zmaxhy attribute values
/*!
\param xminhyParam - the xminhy attribute address - float pointer type
\param xmaxhyParam - the xmaxhy attribute address - float pointer type
\param yminhyParam - the yminhy attribute address - float pointer type
\param ymaxhyParam - the ymaxhy attribute address - float pointer type
\param zminhyParam - the zminhy attribute address - float pointer type
\param zmaxhyParam - the zmaxhy attribute address - float pointer type
*/
void sethydynbound_(float* xminhyParam, float* xmaxhyParam, float* yminhyParam, float* ymaxhyParam, float* zminhyParam, float* zmaxhyParam) {
  hydyn->setBound(*xminhyParam,*xmaxhyParam,*yminhyParam,*ymaxhyParam,*zminhyParam,*zmaxhyParam);
}

//! get ioeos, iozerof attribute values
/*!
\param ioeosParam - the ioeos attribute address - int pointer type
\param iozerofParam - the iozerof attribute address - int pointer type
*/
void gethydynio_(int* ioeosParam, int* iozerofParam) {
  hydyn->getIo(ioeosParam,iozerofParam);
}

//! set ioeos, iozerof attribute values
/*!
\param ioeosParam - the ioeos attribute address - int pointer type
\param iozerofParam - the iozerof attribute address - int pointer type
*/
void sethydynio_(int* ioeosParam, int* iozerofParam) {
  hydyn->setIo(*ioeosParam,*iozerofParam);
}

}
