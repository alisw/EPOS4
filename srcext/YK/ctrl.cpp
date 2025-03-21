//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "ctrl.h"
#include "inc.h"
#include <ctime>
#include <cmath>
#include <iomanip>
//#include <TMath.h>
#include "hdo.h"
#include "ickw.h"
#include "conv.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoCEM.h"
#include "eoBEST.h"
#include "eoPNJL.h"
#include "trancoeff.h"
#include <unistd.h>

#include <execinfo.h>
#include <signal.h>

/*
  EoS functions

     void initXXXX - initialization of the chosen EoS tables with the object „eos”,
                    reading all the variables into the c++ arrays

     void geteosoriginal - the general function reading the variables directly from
                           the tables ε, n B (T, μ B ) , the original ones,
                            if the table is not available, it is created

     void eoshlle - general function, returns all the EoS variables for given ε, n B
*/

using namespace std ;
ofstream ofile ; // cout + cerr output

EoS *eos ;

TransportCoeff *trcoeff ;
Fluid *fluidLarge, *fluidSmall ;
IC_KW *ic_kw ;
Hydro* h ;
time_t start_time, end_time ;

double tau0, tau_max, dtau ;
int maxstep, ienvar ;
//string directory, outkw ;

void actienvar_(int* ienvar_) {ienvar = *ienvar_;}

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  cout<<"Error: signal "<<sig<<endl ;
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

void redirecterrorshlle_(char* filename)
{
    signal(SIGSEGV, handler); // handler for segfault
    ofile.open(filename) ;
    cout.rdbuf(ofile.rdbuf()) ;
    cerr.rdbuf(ofile.rdbuf()) ;
}

void initeoshlle_(char *filename, int* ncols)
{
   // eos = new EoSs(filename,*ncols) ; // << CE ; set BAG_ASYMPT !!
}


void initeoshlle3f_(char *filename, double *B, double *volex0, double *delta0, double *aaa, double *bbb)
{
 eos = new EoS3f(filename, *B, *volex0, *delta0, *aaa, *bbb) ;
// cout << "!!!! EoS substitution: p=e/3\n" ;
// eos = new EoSs("",3);
}

void initeoschiralhlle_(char *fileSmall, char *fileLarge)
{
 eos = new EoSChiral(fileSmall, fileLarge);
}

void initeoshlle1f_(char *filename)
{
    eos = new EoS1f(filename) ;
}

void initbest_(char *filename)
{
    eos = new eoBEST(filename);
}

void initcem_(char* fileMini, char* fileSmall, char* fileLarge)
{
    eos = new eoCEM(fileMini, fileSmall, fileLarge);
}
void initpnjl_(char *filename)
{
    eos = new eoPNJL(filename);
}

void inittrcoeff_(double *etaS, double *zetaS)
{
// trcoeff = new TransportCoeff(0., 0., eos) ;
 trcoeff = new TransportCoeff(*etaS, *zetaS, eos) ;

}

void geteosoriginal_(double *T, double *mu_b, double *mu_q, double *mu_s, double *e, double *n_b, double *n_q, double *n_s, double *p)
{
	eos->eosorginal(*T, *mu_b, *mu_q, *mu_s, *e, *n_b, *n_q, *n_s, *p);
}

void eoshlle_(double *e, double *nb, double *nq, double *ns, double *T, double *mub, double *muq, double *mus, double *p)
{
    eos->eos(*e, *nb, *nq, *ns, *T, *mub, *muq, *mus, *p) ;
}


void eosrangeshlle_(double *emax, double *e0, double *nmax, double *n0, int *ne, int *nn)
{
	eos->eosranges(*emax, *e0, *nmax, *n0, *ne, *nn) ;
}


void initfluidhlle_(int* nx, int* ny, int* nz, double* minx, double* maxx, double* miny, double* maxy, double* minz, double* maxz, double* dtau)
{
    fluidLarge = new Fluid(eos, trcoeff, *nx, *ny, *nz, *minx, *maxx, *miny, *maxy, *minz, *maxz, *dtau) ;
}


void reducegridhlle_(double* dtau, int* factor)
{
 Fluid* f = fluidLarge;
 int n = *factor;
 int nxNew = (f->getNX()+1)/n;
 int nyNew = (f->getNY()+1)/n;
 int nzNew = f->getNZ();
 if(fluidSmall) delete fluidSmall;
 fluidSmall = new Fluid(eos, trcoeff, nxNew, nyNew, nzNew,
   f->getX(0), f->getX(0) + n*f->getDx()*(nxNew-1),
   f->getY(0), f->getY(0) + n*f->getDy()*(nyNew-1),
   f->getZ(0), f->getZ(0) + f->getDz()*(nzNew-1), *dtau);
 for(int ix=0; ix<fluidSmall->getNX(); ix++)
 for(int iy=0; iy<fluidSmall->getNY(); iy++)
 for(int iz=0; iz<fluidSmall->getNZ(); iz++){
  *(fluidSmall->getCell(ix, iy, iz)) = *(f->getCell(n*ix, n*iy, iz));
 }
 for(int ix=0; ix<fluidSmall->getNX(); ix++)
	for(int iy=0; iy<fluidSmall->getNY(); iy++)
	for(int iz=0; iz<fluidSmall->getNZ(); iz++)
	{
  Cell* c = fluidSmall->getCell(ix,iy,iz);
		c->setPrev(X_, fluidSmall->getCell(ix-1, iy, iz)) ;
		c->setNext(X_, fluidSmall->getCell(ix+1, iy, iz)) ;
		c->setPrev(Y_, fluidSmall->getCell(ix, iy-1, iz)) ;
		c->setNext(Y_, fluidSmall->getCell(ix, iy+1, iz)) ;
		c->setPrev(Z_, fluidSmall->getCell(ix, iy, iz-1)) ;
		c->setNext(Z_, fluidSmall->getCell(ix, iy, iz+1)) ;
		c->setPos(ix, iy, iz) ;
	}
 h->setFluid(fluidSmall);
}


void initichlle_(char *filename, double *tau0)
{
    ic_kw = new IC_KW(filename) ; //---------for Klaus
//    ic_kw->setIC(f, eos, *tau0) ; //!!! z-symmetry is included in ic_kw.cpp
}


void icgethlle3f_(double* x, double* y, double* eta, double* e, double* nB, double* nQ, double* nS, double* vx, double* vy, double* vz)
{
	double nu, nd, ns ;
    ic_kw->getICs(*x, *y, *eta, *e, *vx, *vy, *vz, nu, nd, ns) ;
//	if(*e<1.) *e = 0. ;
	*nB = 1./3.*(nu + nd + ns) ;
	*nQ = 1./3.*(2.*nu - nd - ns) ;
	*nS = - ns ;
}


void icgethlle_(double* x, double* y, double* eta, double* e, double* nB, double* nQ, double* nS, double* vx, double* vy, double* vz)
{
	double nu, nd, ns ;
    ic_kw->getICs(*x, *y, *eta, *e, *vx, *vy, *vz, nu, nd, ns) ;
	*nB = *nQ = *nS = 0. ;
}


void icsethlle_(int* ix, int* iy, int* iz, double* tau0, double* e, double* nb, double* nq, double* ns, double* vx, double* vy, double* vz)
{
Cell *c = fluidLarge->getCell(*ix-1,*iy-1,*iz-1) ;
if((*vx)*(*vx)+(*vy)*(*vy)+(*vz)*(*vz)>1.){
//	cerr << "setIC : " << ix << "  " << iy << "  " << iz << "  e = " << e << "  v^2 = " << (*vx)*(*vx)+(*vy)*(*vy)+vz*vz << endl ;
	double factor = sqrt((*vx)*(*vx)+(*vy)*(*vy)+(*vz)*(*vz)) ;
	(*vx) = (*vx)*0.99/factor ;
	(*vy) = (*vy)*0.99/factor ;
	(*vz) = (*vz)*0.99/factor ;
  }
  if(std::isinf(*e) || std::isnan(*e) || std::isinf(*nb) || std::isnan(*nb) || std::isinf(*nq) || std::isnan(*nq) || std::isinf(*ns) || std::isnan(*ns) ||
  std::isinf(*vx) || std::isnan(*vx) || std::isinf(*vy) || std::isnan(*vy) || std::isinf(*vz) || std::isnan(*vz))
  {cout<<"fluid: bad init,"<<setw(14)<<"e"<<setw(14)<<"nb"<<setw(14)<<"nq"<<setw(14)<<"ns"<<setw(14)<<"vx"<<setw(14)<<"vy"<<setw(14)<<"vz"<<endl ;
   cout<<"----------------"<<setw(14)<<e<<setw(14)<<nb<<setw(14)<<nq<<setw(14)<<ns<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<endl ;
   cout<<"cell  "<<ix<<"  "<<iy<<"  "<<iz<<endl ;
   exit(1) ; }
  c->setPrimVar(eos, *tau0, *e, *nb, *nq, *ns, (*vx), (*vy), (*vz)) ;
  c->saveQprev() ;
  if(*e>0.) c->setAllM(1.) ;
}


void inithydrohlle_(double* _tau0, double* _tau_max, double* _dtau)
{
    tau0 = *_tau0 ;
    tau_max = *_tau_max ;
    dtau = *_dtau ;
    h = new Hydro(fluidLarge, eos, trcoeff, tau0, dtau) ;
    maxstep = ceil((tau_max-tau0)/dtau) ;

    start_time = 0;
    time(&start_time);
    h->setNSvalues() ; // initialize viscous terms
    h->setQfull() ; // set Qfull in each cell, in order to output IC correctly
}


void dtauhlle_(double* dtau)
{
  h->setDtau(*dtau) ;
  fluidLarge->setDtau(*dtau);
}


void initoutputhlle_(char* dir)
{
  fluidLarge->initOutput(dir, maxstep, tau0, 2) ;
}


int getmaxstephlle_(void)
{
    return maxstep ;
}


void makestephlle_(int *i)
{
 h->performStep() ;
 if(fluidSmall){  // copy+interpolate fields from small to large grid
  for(int ix=0; ix<fluidLarge->getNX(); ix++)
  for(int iy=0; iy<fluidLarge->getNY(); iy++)
  for(int iz=0; iz<fluidLarge->getNZ(); iz++){
   int ixS = (int)((fluidLarge->getX(ix) - fluidSmall->getX(0))/fluidSmall->getDx());
   int iyS = (int)((fluidLarge->getY(iy) - fluidSmall->getY(0))/fluidSmall->getDy());
   double xm = fluidLarge->getX(ix) - ixS*fluidSmall->getDx() - fluidSmall->getX(0);
   double ym = fluidLarge->getY(iy) - iyS*fluidSmall->getDy() - fluidSmall->getY(0);
   double wx [2] = {1.0 - xm/fluidSmall->getDx(), xm/fluidSmall->getDx()} ;
   double wy [2] = {1.0 - ym/fluidSmall->getDy(), ym/fluidSmall->getDy()} ;
   double Q[7], Qint[7] = {0.,0.,0.,0.,0.,0.,0.};
   double pi_int[4][4], Pi_int=0.;
   for(int i=0; i<4; i++)
   for(int j=0; j<4; j++)
     pi_int[i][j] = 0.0;
   // interpolation loop
   for(int iwx=0; iwx<2; iwx++)
   for(int iwy=0; iwy<2; iwy++){
    Cell *cs = fluidSmall->getCell(ixS+iwx,iyS+iwy,iz);
    cs->getQ(Q);
    for(int i=0; i<7; i++)
     Qint[i] += wx[iwx]*wy[iwy]*Q[i];
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++){
     pi_int[i][j] += wx[iwx]*wy[iwy]*cs->getpi(i,j);
    }
    Pi_int += wx[iwx]*wy[iwy]*cs->getPi();
   } // interpolation loop
   Cell *cl = fluidLarge->getCell(ix,iy,iz);
   cl->saveQprev();
   cl->setQ(Qint);
   for(int i=0; i<4; i++)
   for(int j=0; j<4; j++)
    cl->setpi(i,j,pi_int[i][j]);
   cl->setPi(Pi_int);
  } // end large fluid cell loop
 } // end copy
 if(ienvar==1)fluidLarge->calcEnergies(h->getTau()) ;
}


void getvalueshlle_(int* ix, int* iy, int* iz, double* e, double *p, double *nb, double *nq, double *ns, double* vx, double* vy, double* vz, double* viscCorrCutFlag)
{
    Cell *c = fluidLarge->getCell(*ix-1, *iy-1, *iz-1) ;
    c->getPrimVar(eos, h->getTau(), *e, *p, *nb, *nq, *ns, *vx, *vy, *vz) ;
    *viscCorrCutFlag = c->getViscCorrCutFlag() ;
}


void getvflaghlle_(int* ix, int* iy, int* iz, double* viscCorrCutFlag)
{
  *viscCorrCutFlag = fluidLarge->getCell(*ix-1, *iy-1, *iz-1)->getViscCorrCutFlag() ;
}


void getvischlle_(int* ix, int* iy, int* iz, double *pi, double *Pi)
{
 Cell* c=fluidLarge->getCell(*ix-1, *iy-1, *iz-1) ;
 // fortran and C have reverse array alignment, (a,b) --> [b,a]
 // but since pi[mu][nu] is symmetric, this is not important
 for(int i=0; i<4; i++)
 for(int j=0; j<4; j++){
  *(pi+4*i+j) = c->getpi(i,j) ;
 }
 *Pi=c->getPi();
 //if(c->Pi!=0.0) cout <<"nonzero Pi: " << c->Pi << endl ;
}


void getvisc10hlle_(int* ix, int* iy, int* iz, double* pi, double* Pi)
{
 Cell* c=fluidLarge->getCell(*ix-1, *iy-1, *iz-1) ;
 const double* pix = c->getpi10() ;
 for(int i=0; i<10; i++) pi[i] = pix[i];
 *Pi = c->getPi();
}


double getxhlle_(int *ix)
{
    return fluidLarge->getX(*ix-1) ;
}


double getyhlle_(int *iy)
{
    return fluidLarge->getY(*iy-1) ;
}


double getzhlle_(int *iz)
{
    return fluidLarge->getZ(*iz-1) ;
}


void destroyhlle_(void)
{
//    delete ic ; // do we need IC instance?
    delete fluidLarge ;
    fluidLarge = 0 ;
    if(fluidSmall){
     delete fluidSmall ;
     fluidSmall = 0 ;
    }
    delete h ;
    h = 0 ;
    ofile.close() ; // close output and error file
}


void destroyeoshlle_(void)
{
	delete eos ;
}


double gettimehlle_(void)
{
    end_time=0 ;
    time(&end_time); float diff = difftime(end_time, start_time);
    return diff ;
}

double getenergyhlle_(void)
{
    double ene = 0., e, p, nb, nq, ns, vx, vy, vz, _vx, _vy, _vz ;
    double tau = h->getTau() ;
    for(int ix=0; ix<fluidLarge->getNX(); ix++)
    for(int iy=0; iy<fluidLarge->getNY(); iy++)
    for(int iz=0; iz<fluidLarge->getNZ(); iz++){
	fluidLarge->getCell(ix, iy, iz)->getPrimVar(eos, tau, e, p, nb, nq, ns, _vx, _vy, _vz) ;
	double eta = fluidLarge->getZ(iz) ;
	vx = _vx/(cosh(eta) + _vz*sinh(eta)) ;
	vy = _vy/(cosh(eta) + _vz*sinh(eta)) ;
	vz = (_vz*cosh(eta)+sinh(eta))/(cosh(eta) + _vz*sinh(eta)) ;
	ene += (e+p)/(1.-vx*vx-vy*vy-vz*vz) - p ;
	}
    return e ;
}


void envarget_(double * e, double * efull, double * etotsurf, double * nbsurf)
{
  if(ienvar==1)fluidLarge->getEnergies(*e, *efull, *etotsurf, *nbsurf) ;
  else {e=0;efull=0;etotsurf=0;nbsurf=0;} //ienvar!=1
}
