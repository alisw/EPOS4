#include "Ico.hpp"

//! default class constructor 
Ico::Ico() {}

//! default class destructor 
Ico::~Ico() {
}

// map the fortran common block cchkengy2
extern "C" {
extern struct {
  float esollxx;
  float eistxx;
} cchkengy2_;
}

// get value from fortran common block
void Ico::getFortranEsollxx() { this->esollxx = cchkengy2_.esollxx; }

void Ico::getFortranEistxx() { this->eistxx = cchkengy2_.eistxx; }

// set value from fortran common block
void Ico::setFortranEsollxx() { cchkengy2_.esollxx = this->esollxx; }

void Ico::setFortranEistxx() { cchkengy2_.eistxx = this->eistxx; }

// c++ attribute accessor methods
//! get xmin attribute value
/*!
\return The xmin attribute value - float type
*/
float Ico::getXmin() { return this->xmin; }

//! get xmax attribute value
/*!
\return The xmax attribute value - float type
*/
float Ico::getXmax() { return this->xmax; }

//! get ymin attribute value
/*!
\return The ymin attribute value - float type
*/
float Ico::getYmin() { return this->ymin; }

//! get ymax attribute value
/*!
\return The ymax attribute value - float type
*/
float Ico::getYmax() { return this->ymax; }

//! get zmin attribute value
/*!
\return The zmin attribute value - float type
*/
float Ico::getZmin() { return this->zmin; }

//! get zmax attribute value
/*!
\return The zmax attribute value - float type
*/
float Ico::getZmax() { return this->zmax; }

//! get nx attribute value
/*!
\return The nx attribute value - int type
*/
int Ico::getNx() { return this->nx; }

//! get ny attribute value
/*!
\return The ny attribute value - int type
*/
int Ico::getNy() { return this->ny; }

//! get nz attribute value
/*!
\return The nz attribute value - int type
*/
int Ico::getNz() { return this->nz; }

//! get ee1ico attribute value
/*!
\return The ee1ico attribute value - float type
*/
float Ico::getEe1ico() { return this->ee1ico; }

//! get eistico attribute value
/*!
\return The eistico attribute value - float type
*/
float Ico::getEistico() { return this->eistico; }

//! get ee1hll attribute value
/*!
\return The ee1hll attribute value - float type
*/
float Ico::getEe1hll() { return this->ee1hll; }

//! get ichkengy attribute value
/*!
\return The ichkengy attribute value - int type
*/
int Ico::getIchkengy() { return this->ichkengy; }

//! get esollxx attribute value
/*!
\return The esollxx attribute value - float type
*/
float Ico::getEsollxx() {
  // get value from fortran common block
  this->getFortranEsollxx();
  return this->esollxx;
}

//! get eistxx attribute value
/*!
\return The eistxx attribute value - float type
*/
float Ico::getEistxx() {
  // get value from fortran common block
  this->getFortranEistxx();
  return this->eistxx;
}

//! get energyDensity attribute value
/*!
\return The energyDensity attribute value - Array<float> * type
*/
Array<float> *Ico::getEnergyDensity() { return this->energyDensity; }

//! get flavorDensity attribute value
/*!
\return The flavorDensity attribute value - Array<float> * type
*/
Array<float> *Ico::getFlavorDensity() { return this->flavorDensity; }

//! get velocity attribute value
/*!
\return The velocity attribute value - Array<float> * type
*/
Array<float> *Ico::getVelocity() { return this->velocity; }

//! get energyMomentumTensor attribute value
/*!
\return The energyMomentumTensor attribute value - Array<double> * type
*/
Array<double> *Ico::getEnergyMomentumTensor() { return this->energyMomentumTensor; }

//! get flavorCurrentVector attribute value
/*!
\return The flavorCurrentVector attribute value - Array<double> * type
*/
Array<double> *Ico::getFlavorCurrentVector() { return this->flavorCurrentVector; }

// c++ attribute mutator methods
//! set xmin attribute value
/*!
\param xminParam - the xmin attribute value - float type
*/
void Ico::setXmin(float xminParam) { this->xmin = xminParam; }

//! set xmax attribute value
/*!
\param xmaxParam - the xmax attribute value - float type
*/
void Ico::setXmax(float xmaxParam) { this->xmax = xmaxParam; }

//! set ymin attribute value
/*!
\param yminParam - the ymin attribute value - float type
*/
void Ico::setYmin(float yminParam) { this->ymin = yminParam; }

//! set ymax attribute value
/*!
\param ymaxParam - the ymax attribute value - float type
*/
void Ico::setYmax(float ymaxParam) { this->ymax = ymaxParam; }

//! set zmin attribute value
/*!
\param zminParam - the zmin attribute value - float type
*/
void Ico::setZmin(float zminParam) { this->zmin = zminParam; }

//! set zmax attribute value
/*!
\param zmaxParam - the zmax attribute value - float type
*/
void Ico::setZmax(float zmaxParam) { this->zmax = zmaxParam; }

//! set nx attribute value
/*!
\param nxParam - the nx attribute value - int type
*/
void Ico::setNx(int nxParam) { this->nx = nxParam; }

//! set ny attribute value
/*!
\param nyParam - the ny attribute value - int type
*/
void Ico::setNy(int nyParam) { this->ny = nyParam; }

//! set nz attribute value
/*!
\param nzParam - the nz attribute value - int type
*/
void Ico::setNz(int nzParam) { this->nz = nzParam; }

//! set ee1ico attribute value
/*!
\param ee1icoParam - the ee1ico attribute value - float type
*/
void Ico::setEe1ico(float ee1icoParam) { this->ee1ico = ee1icoParam; }

//! set eistico attribute value
/*!
\param eisticoParam - the eistico attribute value - float type
*/
void Ico::setEistico(float eisticoParam) { this->eistico = eisticoParam; }

//! set ee1hll attribute value
/*!
\param ee1hllParam - the ee1hll attribute value - float type
*/
void Ico::setEe1hll(float ee1hllParam) { this->ee1hll = ee1hllParam; }

//! set ichkengy attribute value
/*!
\param ichkengyParam - the ichkengy attribute value - int type
*/
void Ico::setIchkengy(int ichkengyParam) { this->ichkengy = ichkengyParam; }

//! set esollxx attribute value
/*!
\param esollxxParam - the esollxx attribute value - float type
*/
void Ico::setEsollxx(float esollxxParam) {
  this->esollxx = esollxxParam;
  // set value from fortran common block
  this->setFortranEsollxx();
}

//! set eistxx attribute value
/*!
\param eistxxParam - the eistxx attribute value - float type
*/
void Ico::setEistxx(float eistxxParam) {
  this->eistxx = eistxxParam;
  // set value from fortran common block
  this->setFortranEistxx();
}

//! set energyDensity attribute value
/*!
\param energyDensityParam - the energyDensity attribute value - Array<float> * type
*/
void Ico::setEnergyDensity(Array<float> *energyDensityParam) { this->energyDensity = energyDensityParam; }

//! set flavorDensity attribute value
/*!
\param flavorDensityParam - the flavorDensity attribute value - Array<float> * type
*/
void Ico::setFlavorDensity(Array<float> *flavorDensityParam) { this->flavorDensity = flavorDensityParam; }

//! set velocity attribute value
/*!
\param velocityParam - the velocity attribute value - Array<float> * type
*/
void Ico::setVelocity(Array<float> *velocityParam) { this->velocity = velocityParam; }

//! set energyMomentumTensor attribute value
/*!
\param energyMomentumTensorParam - the energyMomentumTensor attribute value - Array<double> * type
*/
void Ico::setEnergyMomentumTensor(Array<double> *energyMomentumTensorParam) { this->energyMomentumTensor = energyMomentumTensorParam; }

//! set flavorCurrentVector attribute value
/*!
\param flavorCurrentVectorParam - the flavorCurrentVector attribute value - Array<double> * type
*/
void Ico::setFlavorCurrentVector(Array<double> *flavorCurrentVectorParam) { this->flavorCurrentVector = flavorCurrentVectorParam; }

// c++ array attribute methods
//! create energyDensity attribute array of float type
/*!
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/
void Ico::createEnergyDensity(int x, int y, int z) {
  int dim = 3;
  int shape[3] = {x, y, z};
  this->energyDensity = new Array<float>(dim, shape);
}

//! get energyDensity attribute array value
/*!
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\return The energyDensity attribute value - float type
*/
float Ico::getEnergyDensityValue(int x, int y, int z){
  int index[3] = {x, y, z};
  return this->energyDensity->get(index);
}

//! set energyDensity attribute array value
/*!
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\param value - the value to be set - float type
*/
void Ico::setEnergyDensityValue(int x, int y, int z, float value) {
  int index[3] = {x, y, z};
  this->energyDensity->set(index, value);
}

//! create energyDensity initialize method
void Ico::initializeEnergyDensity(float value) {
  this->energyDensity->initialize(value);
}

//! create flavorDensity attribute array of float type
/*!
\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/
void Ico::createFlavorDensity(int n, int x, int y, int z) {
  int dim = 4;
  int shape[4] = {n, x, y, z};
  this->flavorDensity = new Array<float>(dim, shape);
}

//! get flavorDensity attribute array value
/*!
\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\return The flavorDensity attribute value - float type
*/
float Ico::getFlavorDensityValue(int n, int x, int y, int z){
  int index[4] = {n, x, y, z};
  return this->flavorDensity->get(index);
}

//! set flavorDensity attribute array value
/*!
\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\param value - the value to be set - float type
*/
void Ico::setFlavorDensityValue(int n, int x, int y, int z, float value) {
  int index[4] = {n, x, y, z};
  this->flavorDensity->set(index, value);
}

//! create flavorDensity initialize method
void Ico::initializeFlavorDensity(float value) {
  this->flavorDensity->initialize(value);
}

//! create velocity attribute array of float type
/*!
\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/
void Ico::createVelocity(int n, int x, int y, int z) {
  int dim = 4;
  int shape[4] = {n, x, y, z};
  this->velocity = new Array<float>(dim, shape);
}

//! get velocity attribute array value
/*!
\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\return The velocity attribute value - float type
*/
float Ico::getVelocityValue(int n, int x, int y, int z){
  int index[4] = {n, x, y, z};
  return this->velocity->get(index);
}

//! set velocity attribute array value
/*!
\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\param value - the value to be set - float type
*/
void Ico::setVelocityValue(int n, int x, int y, int z, float value) {
  int index[4] = {n, x, y, z};
  this->velocity->set(index, value);
}

//! create velocity initialize method
void Ico::initializeVelocity(float value) {
  this->velocity->initialize(value);
}

//! create energyMomentumTensor attribute array of double type
/*!
\param ip1 - indice - int type
\param ip2 - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/
void Ico::createEnergyMomentumTensor(int ip1, int ip2, int x, int y, int z) {
  int dim = 5;
  int shape[5] = {ip1, ip2, x, y, z};
  this->energyMomentumTensor = new Array<double>(dim, shape);
}

//! get energyMomentumTensor attribute array value
/*!
\param ip1 - indice - int type
\param ip2 - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\return The energyMomentumTensor attribute value - double type
*/
double Ico::getEnergyMomentumTensorValue(int ip1, int ip2, int x, int y, int z){
  int index[5] = {ip1, ip2, x, y, z};
  return this->energyMomentumTensor->get(index);
}

//! set energyMomentumTensor attribute array value
/*!
\param ip1 - indice - int type
\param ip2 - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\param value - the value to be set - double type
*/
void Ico::setEnergyMomentumTensorValue(int ip1, int ip2, int x, int y, int z, double value) {
  int index[5] = {ip1, ip2, x, y, z};
  this->energyMomentumTensor->set(index, value);
}

//! create energyMomentumTensor initialize method
void Ico::initializeEnergyMomentumTensor(double value) {
  this->energyMomentumTensor->initialize(value);
}

//! create flavorCurrentVector attribute array of double type
/*!
\param f - indice - int type
\param ip - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/
void Ico::createFlavorCurrentVector(int f, int ip, int x, int y, int z) {
  int dim = 5;
  int shape[5] = {f, ip, x, y, z};
  this->flavorCurrentVector = new Array<double>(dim, shape);
}

//! get flavorCurrentVector attribute array value
/*!
\param f - indice - int type
\param ip - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\return The flavorCurrentVector attribute value - double type
*/
double Ico::getFlavorCurrentVectorValue(int f, int ip, int x, int y, int z){
  int index[5] = {f, ip, x, y, z};
  return this->flavorCurrentVector->get(index);
}

//! set flavorCurrentVector attribute array value
/*!
\param f - indice - int type
\param ip - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
\param value - the value to be set - double type
*/
void Ico::setFlavorCurrentVectorValue(int f, int ip, int x, int y, int z, double value) {
  int index[5] = {f, ip, x, y, z};
  this->flavorCurrentVector->set(index, value);
}

//! create flavorCurrentVector initialize method
void Ico::initializeFlavorCurrentVector(double value) {
  this->flavorCurrentVector->initialize(value);
}

//! get nx, ny, nz attribute values
/*!
\param nxParam - the nx attribute address - int pointer type
\param nyParam - the ny attribute address - int pointer type
\param nzParam - the nz attribute address - int pointer type
*/
void Ico::getDim(int *nxParam, int *nyParam, int *nzParam) {
  *nxParam = this->getNx();
  *nyParam = this->getNy();
  *nzParam = this->getNz();
}

//! set nx, ny, nz attribute values
/*!
\param nxParam - the nx attribute address - int type
\param nyParam - the ny attribute address - int type
\param nzParam - the nz attribute address - int type
*/
void Ico::setDim(int nxParam, int nyParam, int nzParam) {
  this->setNx(nxParam);
  this->setNy(nyParam);
  this->setNz(nzParam);
}

//! get xmin, xmax, ymin, ymax, zmin, zmax attribute values
/*!
\param xminParam - the xmin attribute address - float pointer type
\param xmaxParam - the xmax attribute address - float pointer type
\param yminParam - the ymin attribute address - float pointer type
\param ymaxParam - the ymax attribute address - float pointer type
\param zminParam - the zmin attribute address - float pointer type
\param zmaxParam - the zmax attribute address - float pointer type
*/
void Ico::getBound(float *xminParam, float *xmaxParam, float *yminParam, float *ymaxParam, float *zminParam, float *zmaxParam) {
  *xminParam = this->getXmin();
  *xmaxParam = this->getXmax();
  *yminParam = this->getYmin();
  *ymaxParam = this->getYmax();
  *zminParam = this->getZmin();
  *zmaxParam = this->getZmax();
}

//! set xmin, xmax, ymin, ymax, zmin, zmax attribute values
/*!
\param xminParam - the xmin attribute address - float type
\param xmaxParam - the xmax attribute address - float type
\param yminParam - the ymin attribute address - float type
\param ymaxParam - the ymax attribute address - float type
\param zminParam - the zmin attribute address - float type
\param zmaxParam - the zmax attribute address - float type
*/
void Ico::setBound(float xminParam, float xmaxParam, float yminParam, float ymaxParam, float zminParam, float zmaxParam) {
  this->setXmin(xminParam);
  this->setXmax(xmaxParam);
  this->setYmin(yminParam);
  this->setYmax(ymaxParam);
  this->setZmin(zminParam);
  this->setZmax(zmaxParam);
}

// Ico instance pointer
Ico *ico;

// Fortran / C wrapper functions
extern "C" {
//! create Ico instance
void createico_() { ico = new Ico(); }

//! delete Ico instance
void destroyico_() { delete ico; }

//! get xmin attribute value
/*!
\return The xmin attribute value - float type
*/
float geticoxmin_() { return ico->getXmin(); }

//! set xmin attribute value
/*!
\param xminParam - the xmin attribute value - float type
*/
void seticoxmin_(float *xminParam) { ico->setXmin(*xminParam); }

//! get xmax attribute value
/*!
\return The xmax attribute value - float type
*/
float geticoxmax_() { return ico->getXmax(); }

//! set xmax attribute value
/*!
\param xmaxParam - the xmax attribute value - float type
*/
void seticoxmax_(float *xmaxParam) { ico->setXmax(*xmaxParam); }

//! get ymin attribute value
/*!
\return The ymin attribute value - float type
*/
float geticoymin_() { return ico->getYmin(); }

//! set ymin attribute value
/*!
\param yminParam - the ymin attribute value - float type
*/
void seticoymin_(float *yminParam) { ico->setYmin(*yminParam); }

//! get ymax attribute value
/*!
\return The ymax attribute value - float type
*/
float geticoymax_() { return ico->getYmax(); }

//! set ymax attribute value
/*!
\param ymaxParam - the ymax attribute value - float type
*/
void seticoymax_(float *ymaxParam) { ico->setYmax(*ymaxParam); }

//! get zmin attribute value
/*!
\return The zmin attribute value - float type
*/
float geticozmin_() { return ico->getZmin(); }

//! set zmin attribute value
/*!
\param zminParam - the zmin attribute value - float type
*/
void seticozmin_(float *zminParam) { ico->setZmin(*zminParam); }

//! get zmax attribute value
/*!
\return The zmax attribute value - float type
*/
float geticozmax_() { return ico->getZmax(); }

//! set zmax attribute value
/*!
\param zmaxParam - the zmax attribute value - float type
*/
void seticozmax_(float *zmaxParam) { ico->setZmax(*zmaxParam); }

//! get nx attribute value
/*!
\return The nx attribute value - int type
*/
int geticonx_() { return ico->getNx(); }

//! set nx attribute value
/*!
\param nxParam - the nx attribute value - int type
*/
void seticonx_(int *nxParam) { ico->setNx(*nxParam); }

//! get ny attribute value
/*!
\return The ny attribute value - int type
*/
int geticony_() { return ico->getNy(); }

//! set ny attribute value
/*!
\param nyParam - the ny attribute value - int type
*/
void seticony_(int *nyParam) { ico->setNy(*nyParam); }

//! get nz attribute value
/*!
\return The nz attribute value - int type
*/
int geticonz_() { return ico->getNz(); }

//! set nz attribute value
/*!
\param nzParam - the nz attribute value - int type
*/
void seticonz_(int *nzParam) { ico->setNz(*nzParam); }

//! get ee1ico attribute value
/*!
\return The ee1ico attribute value - float type
*/
float geticoee1ico_() { return ico->getEe1ico(); }

//! set ee1ico attribute value
/*!
\param ee1icoParam - the ee1ico attribute value - float type
*/
void seticoee1ico_(float *ee1icoParam) { ico->setEe1ico(*ee1icoParam); }

//! get eistico attribute value
/*!
\return The eistico attribute value - float type
*/
float geticoeistico_() { return ico->getEistico(); }

//! set eistico attribute value
/*!
\param eisticoParam - the eistico attribute value - float type
*/
void seticoeistico_(float *eisticoParam) { ico->setEistico(*eisticoParam); }

//! get ee1hll attribute value
/*!
\return The ee1hll attribute value - float type
*/
float geticoee1hll_() { return ico->getEe1hll(); }

//! set ee1hll attribute value
/*!
\param ee1hllParam - the ee1hll attribute value - float type
*/
void seticoee1hll_(float *ee1hllParam) { ico->setEe1hll(*ee1hllParam); }

//! get ichkengy attribute value
/*!
\return The ichkengy attribute value - int type
*/
int geticoichkengy_() { return ico->getIchkengy(); }

//! set ichkengy attribute value
/*!
\param ichkengyParam - the ichkengy attribute value - int type
*/
void seticoichkengy_(int *ichkengyParam) { ico->setIchkengy(*ichkengyParam); }

//! get esollxx attribute value
/*!
\return The esollxx attribute value - float type
*/
float geticoesollxx_() { return ico->getEsollxx(); }

//! set esollxx attribute value
/*!
\param esollxxParam - the esollxx attribute value - float type
*/
void seticoesollxx_(float *esollxxParam) { ico->setEsollxx(*esollxxParam); }

//! get eistxx attribute value
/*!
\return The eistxx attribute value - float type
*/
float geticoeistxx_() { return ico->getEistxx(); }

//! set eistxx attribute value
/*!
\param eistxxParam - the eistxx attribute value - float type
*/
void seticoeistxx_(float *eistxxParam) { ico->setEistxx(*eistxxParam); }

//! create energyDensity array attribute
/*!\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/void createicoicoe_(int* x, int* y, int* z) {
  ico->createEnergyDensity(*x, *y, *z);
}

//! initialize energyDensity array attribute
void initializeicoicoe_(float *value) {
  ico->initializeEnergyDensity(*value);
}

//! set energyDensity array value
/*!
\param x - indice - int pointer type
\param y - indice - int pointer type
\param z - indice - int pointer type
\param value - value to set - float pointer type
*/
void seticoicoe_(int *x, int *y, int *z, float *value) {
  ico->setEnergyDensityValue(*x, *y, *z, *value);
}

//! get energyDensity array value
/*!
\param x indice - int pointer type
\param y indice - int pointer type
\param z indice - int pointer type
\return value - value to get - float pointer type
*/
float geticoicoe_(int *x, int *y, int *z) {
  return ico->getEnergyDensityValue(*x, *y, *z);
}

//! delete energyDensity array
void destroyicoicoe_() { delete ico->getEnergyDensity(); }

//! create flavorDensity array attribute
/*!\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/void createicoicof_(int* n, int* x, int* y, int* z) {
  ico->createFlavorDensity(*n, *x, *y, *z);
}

//! initialize flavorDensity array attribute
void initializeicoicof_(float *value) {
  ico->initializeFlavorDensity(*value);
}

//! set flavorDensity array value
/*!
\param n - indice - int pointer type
\param x - indice - int pointer type
\param y - indice - int pointer type
\param z - indice - int pointer type
\param value - value to set - float pointer type
*/
void seticoicof_(int *n, int *x, int *y, int *z, float *value) {
  ico->setFlavorDensityValue(*n, *x, *y, *z, *value);
}

//! get flavorDensity array value
/*!
\param n indice - int pointer type
\param x indice - int pointer type
\param y indice - int pointer type
\param z indice - int pointer type
\return value - value to get - float pointer type
*/
float geticoicof_(int *n, int *x, int *y, int *z) {
  return ico->getFlavorDensityValue(*n, *x, *y, *z);
}

//! delete flavorDensity array
void destroyicoicof_() { delete ico->getFlavorDensity(); }

//! create velocity array attribute
/*!\param n - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/void createicoicov_(int* n, int* x, int* y, int* z) {
  ico->createVelocity(*n, *x, *y, *z);
}

//! initialize velocity array attribute
void initializeicoicov_(float *value) {
  ico->initializeVelocity(*value);
}

//! set velocity array value
/*!
\param n - indice - int pointer type
\param x - indice - int pointer type
\param y - indice - int pointer type
\param z - indice - int pointer type
\param value - value to set - float pointer type
*/
void seticoicov_(int *n, int *x, int *y, int *z, float *value) {
  ico->setVelocityValue(*n, *x, *y, *z, *value);
}

//! get velocity array value
/*!
\param n indice - int pointer type
\param x indice - int pointer type
\param y indice - int pointer type
\param z indice - int pointer type
\return value - value to get - float pointer type
*/
float geticoicov_(int *n, int *x, int *y, int *z) {
  return ico->getVelocityValue(*n, *x, *y, *z);
}

//! delete velocity array
void destroyicoicov_() { delete ico->getVelocity(); }

//! create energyMomentumTensor array attribute
/*!\param ip1 - indice - int type
\param ip2 - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/void createicoicot_(int* ip1, int* ip2, int* x, int* y, int* z) {
  ico->createEnergyMomentumTensor(*ip1, *ip2, *x, *y, *z);
}

//! initialize energyMomentumTensor array attribute
void initializeicoicot_(double *value) {
  ico->initializeEnergyMomentumTensor(*value);
}

//! set energyMomentumTensor array value
/*!
\param ip1 - indice - int pointer type
\param ip2 - indice - int pointer type
\param x - indice - int pointer type
\param y - indice - int pointer type
\param z - indice - int pointer type
\param value - value to set - double pointer type
*/
void seticoicot_(int *ip1, int *ip2, int *x, int *y, int *z, double *value) {
  ico->setEnergyMomentumTensorValue(*ip1, *ip2, *x, *y, *z, *value);
}

//! get energyMomentumTensor array value
/*!
\param ip1 indice - int pointer type
\param ip2 indice - int pointer type
\param x indice - int pointer type
\param y indice - int pointer type
\param z indice - int pointer type
\return value - value to get - double pointer type
*/
double geticoicot_(int *ip1, int *ip2, int *x, int *y, int *z) {
  return ico->getEnergyMomentumTensorValue(*ip1, *ip2, *x, *y, *z);
}

//! delete energyMomentumTensor array
void destroyicoicot_() { delete ico->getEnergyMomentumTensor(); }

//! create flavorCurrentVector array attribute
/*!\param f - indice - int type
\param ip - indice - int type
\param x - indice - int type
\param y - indice - int type
\param z - indice - int type
*/void createicoicoc_(int* f, int* ip, int* x, int* y, int* z) {
  ico->createFlavorCurrentVector(*f, *ip, *x, *y, *z);
}

//! initialize flavorCurrentVector array attribute
void initializeicoicoc_(double *value) {
  ico->initializeFlavorCurrentVector(*value);
}

//! set flavorCurrentVector array value
/*!
\param f - indice - int pointer type
\param ip - indice - int pointer type
\param x - indice - int pointer type
\param y - indice - int pointer type
\param z - indice - int pointer type
\param value - value to set - double pointer type
*/
void seticoicoc_(int *f, int *ip, int *x, int *y, int *z, double *value) {
  ico->setFlavorCurrentVectorValue(*f, *ip, *x, *y, *z, *value);
}

//! get flavorCurrentVector array value
/*!
\param f indice - int pointer type
\param ip indice - int pointer type
\param x indice - int pointer type
\param y indice - int pointer type
\param z indice - int pointer type
\return value - value to get - double pointer type
*/
double geticoicoc_(int *f, int *ip, int *x, int *y, int *z) {
  return ico->getFlavorCurrentVectorValue(*f, *ip, *x, *y, *z);
}

//! delete flavorCurrentVector array
void destroyicoicoc_() { delete ico->getFlavorCurrentVector(); }

//! get nx, ny, nz attribute values
/*!
\param nxParam - the nx attribute address - int pointer type
\param nyParam - the ny attribute address - int pointer type
\param nzParam - the nz attribute address - int pointer type
*/
void geticodim_(int* nxParam, int* nyParam, int* nzParam) {
  ico->getDim(nxParam,nyParam,nzParam);
}

//! set nx, ny, nz attribute values
/*!
\param nxParam - the nx attribute address - int pointer type
\param nyParam - the ny attribute address - int pointer type
\param nzParam - the nz attribute address - int pointer type
*/
void seticodim_(int* nxParam, int* nyParam, int* nzParam) {
  ico->setDim(*nxParam,*nyParam,*nzParam);
}

//! get xmin, xmax, ymin, ymax, zmin, zmax attribute values
/*!
\param xminParam - the xmin attribute address - float pointer type
\param xmaxParam - the xmax attribute address - float pointer type
\param yminParam - the ymin attribute address - float pointer type
\param ymaxParam - the ymax attribute address - float pointer type
\param zminParam - the zmin attribute address - float pointer type
\param zmaxParam - the zmax attribute address - float pointer type
*/
void geticobound_(float* xminParam, float* xmaxParam, float* yminParam, float* ymaxParam, float* zminParam, float* zmaxParam) {
  ico->getBound(xminParam,xmaxParam,yminParam,ymaxParam,zminParam,zmaxParam);
}

//! set xmin, xmax, ymin, ymax, zmin, zmax attribute values
/*!
\param xminParam - the xmin attribute address - float pointer type
\param xmaxParam - the xmax attribute address - float pointer type
\param yminParam - the ymin attribute address - float pointer type
\param ymaxParam - the ymax attribute address - float pointer type
\param zminParam - the zmin attribute address - float pointer type
\param zmaxParam - the zmax attribute address - float pointer type
*/
void seticobound_(float* xminParam, float* xmaxParam, float* yminParam, float* ymaxParam, float* zminParam, float* zmaxParam) {
  ico->setBound(*xminParam,*xmaxParam,*yminParam,*ymaxParam,*zminParam,*zmaxParam);
}

}
