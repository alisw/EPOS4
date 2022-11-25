//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <TError.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <TF1.h>
#include <time.h>
#include "dbs.h"
#include "uku.h"

#include <iostream>
#include <fstream>

#include "es.h"

using namespace std ;

ofstream m2file ;

DatabasePDG2 *database ;
TC *tc ;

const int cStringLength = 30 ;
double dummy = 0.2 ; //dummy variable for random number generator call

        const double gevtofm = 5.067728853 ;
#define S_OK 0
#define S_NOROOT 1

        int eosChoice = EOS3F ;

        double gammaS = 0.93 ;
        double B = 0.38 ;
        double mS = 0.12 ;

        double volex0 = 0 ;
        double delta0 = 1e-6 ;
        double aaa = 0 ;
        double bbb = 0 ;
        double muc = 0.2 ;

        double emax = 40. ;
        double e0 = 0.01 ;
        double anmax = 0. ;
        double bnmax = 0. ;
        double an0 = 0.0001 ;
        int nxe = 65 ; 
        int nxn = 33 ;//KW


        double xe_max = log(emax/e0+1) ;
        double del_xe =    xe_max/(nxe-1) ;

        double muBmax = 0.4 ;
        double muQmax = 0.13 ;
        double muSmax = 0.35 ;
        double muBmx = 0.4 ;
        double muQmx = 0.13 ;
        double muSmx = 0.35 ;
        double temmax = 1.2 ;
        int maxIter = 100 ;
        double Taccuracy = 0.000005 ;

        double Raccuracy = 1e-4 ;

void checkroot(double (*f)(double, double*), double* par, double xmin, double xmax, int *n);
double findroot(double (*f)(double, double*), double* par, double xmin, double xmax, int &status);
void mix(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p);

void SolveMix2(double e, double nb, double nq, double ns, double &T, double &mub, double &muq, double &mus, double &p) ;
void solve2(double e, double nb, double &T, double &mub, double &p) ;

void seteostype_(int *type)
{
        eosChoice = *type ;
}

extern "C"{
 double dranf_(double*) ;
 void ranfini_(double*, int*, int*) ;
 double drangen_(double*) ;
}

void loadptldky_(char *fileparticles, char *filedecay)
{

        database = new DatabasePDG2(fileparticles, filedecay);
         
        if(!database->LoadData()) exit(1) ;
        
        database->SetMassRange(0.05, 10.0); //-------without PHOTONS
        
        database->SetWidthRange(0., 10.);

        database->SortDecayingResonances() ;
        //database->DumpDecays();

//        for(int i=0; i<database->GetNParticles(); i++){
//                 ParticlePDG2 *particle = database->GetPDGParticleByIndex(i);
//                double mass = particle->GetMass() ;
//                if(particle->GetNDecayChannels()==0.) cout << "0br: " << setw(14) << particle->GetPDG() << setw(24) << particle->GetName() 
//                << setw(14) << particle->GetMass() << setw(14) << particle->GetWidth() << endl ;
//        }
        
        //database->DumpData() ;
}

void destroyptldky_(void)
{
        delete database ;
}


extern "C" void checkdecays_(void)
{
        for(int i=0; i<database->GetNParticles(); i++){
                 ParticlePDG2 *particle = database->GetPDGParticleByIndex(i);
                if(particle->GetNDecayChannels()==0.) cout << "no_decay: " << setw(14) << particle->GetPDG() << setw(24) << particle->GetName() 
                << setw(14) << particle->GetMass() << setw(14) << particle->GetWidth() << endl ;
        }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void setparameters_(double *_gammaS, double *_B, double *_mS, double *_volex0, double *_delta0, double *_aaa, double *_bbb, double *_muc, double *_emax, double *_e0,
        double *_anmax, double *_bnmax, double *_an0, int *_nxe, int *_nxn, double *_muBmax, double *_muQmax, double *_muSmax,   double *_muBmx, double *_muQmx, double *_muSmx, 
        double *_temmax, int *_maxIter, double *_Taccuracy)
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
{
        gammaS = *_gammaS ;
        B = *_B ;
        mS = *_mS ;

        volex0 = *_volex0 ;     
        delta0 = *_delta0 ;
        aaa = *_aaa ;
        bbb = *_bbb ;
        muc = *_muc ;

        emax  = *_emax ;
        anmax = *_anmax ;
        bnmax = *_bnmax ;
        e0 = *_e0 ;
        an0 = *_an0 ;
        xe_max = log(emax/e0+1) ;
        nxe = *_nxe ;
        nxn = *_nxn ;
        del_xe = xe_max/(nxe-1) ;

        muBmax = *_muBmax ;
        muQmax = *_muQmax ;
        muSmax = *_muSmax ;
        muBmx = *_muBmx ;
        muQmx = *_muQmx ;
        muSmx = *_muSmx ;
        temmax = *_temmax ;
        maxIter = *_maxIter ;
        Taccuracy = *_Taccuracy ;
                
        //cout <<"++++++++++setparameters:  aaa = "<<aaa<<"  bbb = "<<bbb<<endl;
        //cout << "+++++++ B  = "<<B<< "   volex0 = "<<volex0<<"   delta0 = "<<delta0 <<endl;
        
}

void eosihlle_(double *T, double *mu_b, double *mu_q, double *mu_s, double *e, double *n_b, double *n_q, double *n_s, double *p) 
{
       mix(*T, *mu_b, *mu_q, *mu_s, *e, *n_b, *n_q, *n_s, *p);
}

void eosohlle_(double *e, double *n_b, double *n_q, double *n_s, double *T, double *mu_b, double *mu_q, double *mu_s, double *p)

{
        if(eosChoice==EOS3F)
        SolveMix2(*e, *n_b, *n_q, *n_s, *T, *mu_b, *mu_q, *mu_s, *p) ;
        else{
        solve2(*e, *n_b, *T, *mu_b, *p) ;
        *mu_q = *mu_s = 0. ;
        }
}

double sgn(double x)
{
        if(x>0) return 1. ;
        else if(x<0.) return -1. ;
        else return 0. ;
}


void HG_boltzmann(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{
        p = 0. ;
        n_b = 0. ;
        n_q = 0. ;
        n_s = 0. ;
        e = 0. ;
        if(T<0.005) return ;
        for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
        ParticlePDG2 *particle = database->GetPDGParticleByIndex(currParticle);
        int pid = particle->GetPDG();
        if(abs(pid)==310 || abs(pid)==130)
         continue;
        Double_t g = 2.*particle->GetSpin() + 1.;                                    // degeneracy factor
        Double_t m = particle->GetMass();                                                // PDG table mass
        double mu = particle->GetBaryonNumber()*mu_b + particle->GetElectricCharge()*mu_q + particle->GetStrangeness()*mu_s ;
        if(mu-m>0.) mu = m ;
        double n = g/(2.*pow(TMath::Pi(),2))*m*m*T*exp(mu/T)
                *TMath::BesselK(2,m/T)*pow(gammaS,particle->GetStrangeQNumber())*pow(gevtofm, 3) ;
        p += n*T ;
//        cout << setw(8) << particle->GetPDG() << setw(10) << mu << setw(15) << exp(mu/T) << setw(15) << n*T << endl ;
        n_b += n*particle->GetBaryonNumber() ;
        n_q += n*particle->GetElectricCharge() ;
        n_s += n*particle->GetStrangeness() ;
        e += g/(2.*pow(TMath::Pi(),2))*m*m*T*(3.*T*TMath::BesselK(2,m/T)+m*TMath::BesselK(1,m/T))
                *exp(mu/T)*pow(gammaS,particle->GetStrangeQNumber())*pow(gevtofm, 3) ;

//        if(particle->GetPDG()==2212) cout << "primary protons : " << n[currParticle] << endl ;
        }
}


// boltzmann approximation + excluded volume correction
void HG_exv2(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{
        double p0 = 0. ;
        double n = 0. ;
        p = 0. ;
        n_b = 0. ;
        n_q = 0. ;
        n_s = 0. ;
        e = 0. ;
        if(T<0.005) return ;

        //KW to avoid inf
        double excs=0.;
        double mubb=mu_b, muqq=mu_q, muss=mu_s;
        if(mubb> muBmx){excs=excs+mubb-muBmx;mubb= muBmx;}
        if(mubb<-muBmx){excs=excs-mubb-muBmx;mubb=-muBmx;}
        if(muqq> muQmx){excs=excs+muqq-muQmx;muqq= muQmx;}
        if(muqq<-muQmx){excs=excs-muqq-muQmx;muqq=-muQmx;}
        if(muss> muSmx){excs=excs+muss-muSmx;muss= muSmx;}
        if(muss<-muSmx){excs=excs-muss-muSmx;muss=-muSmx;}

        double v =  fabs(volex0) ; //1.43676 ; // corresponds to 4/3 pi r^3, r=0.7 fm

        // first we determine p0 - without excluded volume correction
        for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
        ParticlePDG2 *particle = database->GetPDGParticleByIndex(currParticle);
        int pid = particle->GetPDG();
        if(abs(pid)==310 || abs(pid)==130)
         continue;
        Double_t g = 2.*particle->GetSpin() + 1.;                                    // degeneracy factor
        Double_t m = particle->GetMass();                                                // PDG table mass
        double mu = particle->GetBaryonNumber()*mubb + particle->GetElectricCharge()*muqq + particle->GetStrangeness()*muss ;
        if(mu>700.*T) mu = 700.*T;  //new 3118

        double dp = TMath::BesselK(2,m/T)*exp(mu/T) ;

        dp *= g/(2.*pow(TMath::Pi(),2))*m*m*T*T*pow(gammaS,particle->GetStrangeQNumber())*pow(gevtofm, 3) ;
        p0 += dp ;
        }

        //cout<<"KWtest   "<<p0<<"  "<< excs/T <<"  "<< exp(200)<<endl;
        p0 = p0 * exp( fmin(200.,excs/T) );
        //cout<<"KWtest   "<<p0<<endl;

        double pmin = 0., pmax = p0 ;
        do{
                p = (pmin+pmax)/2. ;
                if( p0*exp(-v*p/T) - p > 0 ) pmin = p ; else pmax = p ;
        }while((pmax-pmin)/pmax>0.0001) ;
        p = (pmin+pmax)/2. ;

        // energy density, n's
        for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
        ParticlePDG2 *particle = database->GetPDGParticleByIndex(currParticle);
        int pid = particle->GetPDG();
        if(abs(pid)==310 || abs(pid)==130)
         continue;
        Double_t g = 2.*particle->GetSpin() + 1.;                                    // degeneracy factor
        Double_t m = particle->GetMass();                                                // PDG table mass
        double mu = particle->GetBaryonNumber()*mu_b + particle->GetElectricCharge()*mu_q + particle->GetStrangeness()*mu_s ;
        if(mu>700.*T) mu = 700.*T;  //new 3118

        double dn = TMath::BesselK(2,m/T)*exp(mu/T-v*p/T) ;
        double de = (3.*T*TMath::BesselK(2,m/T)+m*TMath::BesselK(1,m/T))*exp(mu/T-v*p/T) ;

        dn *= g/(2.*pow(TMath::Pi(),2))*m*m*T*pow(gammaS,particle->GetStrangeQNumber())*pow(gevtofm, 3) ;
        de *= g/(2.*pow(TMath::Pi(),2))*m*m*T*pow(gammaS,particle->GetStrangeQNumber())*pow(gevtofm, 3) ;
        n_b += dn*particle->GetBaryonNumber() ;
        n_q += dn*particle->GetElectricCharge() ;
        n_s += dn*particle->GetStrangeness() ;
        e += de ;
        n += dn ;
        }
        if(volex0 > 0)        e = e/(1.+v*n) ;  //volex < 0 for testing
        n_b = n_b/(1.+v*n) ;
        n_q = n_q/(1.+v*n) ;
        n_s = n_s/(1.+v*n) ;
}



double fhadron_p(double *x, double *par)
{
        double &T = par[0] ;
        double &mu = par[1] ;
        double &m = par[2] ;
        double &sign = par[3] ;
        double &k = x[0] ;
        return -sign*k*k*log( 1 - sign * exp( - sqrt(k*k+m*m)/T + mu/T ) ) ;
}


double fhadron_epsilon(double *x, double *par)
{
        double &T = par[0] ;
        double &mu = par[1] ;
        double &m = par[2] ;
        double &sign = par[3] ;
        double &k = x[0] ;
        return k*k*sqrt(k*k+m*m)/(exp( sqrt(k*k+m*m)/T - mu/T ) - sign ) ;
}


double fhadron_n(double *x, double *par)
{
        double &T = par[0] ;
        double &mu = par[1] ;
        double &m = par[2] ;
        double &sign = par[3] ;
        double &k = x[0] ;
        return k*k/(exp( sqrt(k*k+m*m)/T - mu/T ) - sign ) ;
}


void HG(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{
        p = 0. ;
        n_b = 0. ;
        n_q = 0. ;
        n_s = 0. ;
        e = 0. ;
        if(T<0.005) return ;
        for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
        ParticlePDG2 *particle = database->GetPDGParticleByIndex(currParticle);
        int pid = particle->GetPDG();
        if(abs(pid)==310 || abs(pid)==130)
         continue;
        Double_t g = 2.*particle->GetSpin() + 1.;                                    // degeneracy factor
        Double_t m = particle->GetMass();                                                // PDG table mass
        double sign = int(2.*particle->GetSpin()) & 1 ? -1. : 1. ;
        double mu = particle->GetBaryonNumber()*mu_b + particle->GetElectricCharge()*mu_q + particle->GetStrangeness()*mu_s ;
//        if(mu-m>0.) mu = m ;

        double dn=0., dp=0., de=0. ;
        // pressure
        TF1 *fp = new TF1("fp", fhadron_p, 0., 10., 4) ;
        fp->SetParameters(T, mu, m, sign) ;
        dp = g*T/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*fp->Integral(0., 10.) ;

        // energy density
        TF1 *feps = new TF1("feps", fhadron_epsilon, 0., 10., 4) ;
        feps->SetParameters(T, mu, m, sign) ;
        de = g/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*feps->Integral(0., 10.) ;

        // numbers
        TF1 *fn = new TF1("fn", fhadron_n, 0., 10., 4) ;
        fn->SetParameters(T, mu, m, sign) ;
        dn = g/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*fn->Integral(0., 10.) ;

        delete fp ;
        delete feps ;
        delete fn ;

        n_b += dn*particle->GetBaryonNumber() ;
        n_q += dn*particle->GetElectricCharge() ;
        n_s += dn*particle->GetStrangeness() ;
        p += dp ;
        e += de ;
        }
}


void qs_contribution(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) ;


void QGP(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{
        const double muU = mu_b/3. + 2./3.*mu_q ;
        const double muD = mu_b/3. - 1./3.*mu_q ;
        if(T<0.005) { e = n_b = n_q = n_s = p = 0.; return ; }
        const double g_l = 6. ;
        const double g_g = 16 ;

        p = g_l/6./pow(TMath::Pi(),2)*(pow(muU,4)/4. + pow(TMath::Pi()*muU*T,2)/2. + 7.*pow(TMath::Pi()*T,4)/60.) +
                g_l/6./pow(TMath::Pi(),2)*(pow(muD,4)/4. + pow(TMath::Pi()*muD*T,2)/2. + 7.*pow(TMath::Pi()*T,4)/60.) +
                g_g*pow(TMath::Pi()*T*T,2)/90. ;
        n_b = g_l/6./pow(TMath::Pi(),2)*(pow(muU,3) + pow(TMath::Pi(),2)*muU*T*T + pow(muD,3) + pow(TMath::Pi(),2)*muD*T*T)/3. ;
        n_q = g_l/6./pow(TMath::Pi(),2)*(2.*pow(muU,3) + 2.*pow(TMath::Pi(),2)*muU*T*T - pow(muD,3) - pow(TMath::Pi(),2)*muD*T*T)/3. ;
        n_s = 0. ;
        p *= pow(gevtofm, 3) ;
        e = 3.*p + B ;
        p -= B ;
        n_b *= pow(gevtofm, 3) ;
        n_q *= pow(gevtofm, 3) ;
        n_s *= pow(gevtofm, 3) ;
        // s-quark contributions
        double es=0., nbs=0., nqs=0., nss=0., ps=0. ;
        qs_contribution(T, mu_b, mu_q, mu_s, es, nbs, nqs, nss, ps) ; 
        e += es ;
        n_b += nbs ;
        n_q += nqs ;
        n_s += nss ;
        p += ps ;
}


double fqs_p(double *x, double *par)
{
        double &T = par[0] ;
        double &mu = par[1] ;
        double &k = x[0] ;
        return k*k*log( 1 + exp( - sqrt(k*k+mS*mS)/T + mu/T ) ) ;
}


double fqs_epsilon(double *x, double *par)
{
        double &T = par[0] ;
        double &mu = par[1] ;
        double &k = x[0] ;
        return k*k*sqrt(k*k+mS*mS)/(exp( sqrt(k*k+mS*mS)/T - mu/T ) + 1) ;
}


double fqs_n(double *x, double *par)
{
        double &T = par[0] ;
        double &mu = par[1] ;
        double &k = x[0] ;
        return k*k/(exp( sqrt(k*k+mS*mS)/T - mu/T ) + 1) ;
}


void qs_contribution(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{
        const double g_l = 6. ;
        const double mu = mu_b/3. - mu_s - 1./3.*mu_q ;

        // pressure
        TF1 *fp = new TF1("fp", fqs_p, 0., 10., 2) ;
        fp->SetParameters(T, mu) ;
        p = g_l*T/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*fp->Integral(0., 100.) ;
        fp->SetParameters(T, -mu) ;
        p += g_l*T/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*fp->Integral(0., 100.) ;

        // energy density
        TF1 *feps = new TF1("feps", fqs_epsilon, 0., 10., 2) ;
        feps->SetParameters(T, mu) ;
        e = g_l/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*feps->Integral(0., 100.) ;
        feps->SetParameters(T, -mu) ;
        e += g_l/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*feps->Integral(0., 100.) ;


        // numbers
        TF1 *fn = new TF1("fn", fqs_n, 0., 10., 2) ;
        fn->SetParameters(T, mu) ;
        double n = g_l/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*fn->Integral(0., 100.) ;
        n_b = 1./3.*n ;
        n_q = -1./3.*n ;
        n_s = -n ;
        fn->SetParameters(T, -mu) ;
        n = g_l/2./pow(TMath::Pi(),2)*pow(gevtofm, 3)*fn->Integral(0., 100.) ;
        n_b += -1./3.*n ;
        n_q += 1./3.*n ;
        n_s += n ;

        delete fp ;
        delete feps ;
        delete fn ;
}

void QGP_uds(double T, double muU, double muD, double muS, double &e, double& n_u, double& n_d, double& n_s, double& p)
{
        const double g_l = 6. ;
        const double g_g = 16 ;

        double deltap = 0.;
        p = g_l/6./pow(TMath::Pi(),2)*(pow(muU,4)/4. + pow(TMath::Pi()*muU*T,2)/2. + 7.*pow(TMath::Pi()*T,4)/60.) +
                g_l/6./pow(TMath::Pi(),2)*(pow(muS,4)/4. + pow(TMath::Pi()*muS*T,2)/2. + 7.*pow(TMath::Pi()*T,4)/60.) +
                g_g*pow(TMath::Pi()*T*T,2)/90. ;
        n_u = g_l/6./pow(TMath::Pi(),2)*(pow(muU,3) + pow(TMath::Pi(),2)*muU*T*T) ;
        n_d = g_l/6./pow(TMath::Pi(),2)*(pow(muD,3) + pow(TMath::Pi(),2)*muD*T*T) ;
        n_s = 0. ;
        for(int i=1; i<2; i++){
                p += g_l/2./pow(TMath::Pi(),2)*mS*mS*T*T*pow(-1.,i+1)/(i*i)*TMath::BesselK(2, i*mS/T)*(exp(i*muS/T)+exp(-i*muS/T)) ;
                n_s += g_l/2./pow(TMath::Pi(),2)*mS*mS*T*pow(-1.,i+1)/i*TMath::BesselK(2, i*mS/T)*(-exp(i*muS/T)+exp(-i*muS/T)) ;
                deltap += g_l/2./pow(TMath::Pi(),2)*mS*mS*mS*T*pow(-1.,i+1)/i*TMath::BesselK(1, i*mS/T)*(exp(i*muS/T)+exp(-i*muS/T)) ;
        }
        p *= pow(gevtofm, 3) ;
        e = 3.*p + B + deltap*pow(gevtofm, 3) ;
        p -= B ;
        n_u *= pow(gevtofm, 3) ;
        n_d *= pow(gevtofm, 3) ;
        n_s *= pow(gevtofm, 3) ;
}

double getT_c(double mu_b, double mu_q, double mu_s)
{
     double eH, nbH, nqH, nsH, pH, eQ, nbQ, nqQ, nsQ, pQ ;
     double T, Tmin = 0.05, Tmax = 0.18 ;
     do{
     T = 0.5*(Tmin+Tmax) ;
     HG_exv2(T, mu_b, mu_q, mu_s, eH, nbH, nqH, nsH, pH) ;
     QGP(T, mu_b, mu_q, mu_s, eQ, nbQ, nqQ, nsQ, pQ) ;
     if(pQ>pH) Tmax = T ;
     else Tmin = T ;
     }while((Tmax-Tmin)>0.00001) ;
     return 0.5*(Tmin+Tmax) ;
}

void initeos3tc_(char *filename) { tc = new TC(filename) ; }

TC::TC(char *filename){
        ifstream fin (filename) ;
        if(!fin) { cout << "cannot open Tc file; exiting\n" ; exit(0) ; }
        fin >> aa >> bb >> nn ;
        if(fabs(aa-B)     >1e-5){ cout <<"\n\n\n      B =      "<<aa<<" instead of "<<B     <<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(bb-volex0)>1e-5){ cout <<"\n\n\n      volex0 = "<<bb<<" instead of "<<volex0<<"; exiting\n\n\n\n" ; exit(0) ; }
        //cout <<"++++++++++ B = "<<B<<"  delta0 = "<<delta0<<endl;
        TcTab = new double [nn*nn*nn] ;
        for(int ib=0; ib<nn; ib++)
        for(int iq=0; iq<nn; iq++)
        for(int is=0; is<nn; is++) {fin >> TcTab[index(ib,iq,is)] ;}
        fin.close() ;}
TC::~TC(){delete [] TcTab ;}

void TC::getu(double mu, double muMax, int &i, double &u){
        const double dmu = 2*muMax/(nn-1) ;
        i = (int) ((mu+muMax)/dmu) ;
        if(i>nn-2) i = nn - 2 ;
        if(i<0) i = 0 ;
        double ddmu = (mu+muMax) - i*dmu ;
        u = ddmu/dmu ;
}

double TC::pi(double mub, double muq, double mus){
        int ib , is , iq ;
        double ub, uq, us ;
        getu(mub, muBmx, ib, ub);
        getu(muq, muQmx, iq, uq);
        getu(mus, muSmx, is, us);
        double wb [2] = {1.-ub, ub} ;
        double wq [2] = {1.-uq, uq} ;
        double ws [2] = {1.-us, us} ;

        double T = 0. ;
        for(int jb=0; jb<2; jb++)
        for(int jq=0; jq<2; jq++)
        for(int js=0; js<2; js++){T += wb[jb]*wq[jq]*ws[js]*TcTab[index(ib+jb,iq+jq,is+js)] ;}
        // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
        return T;
}

void maketctable_(char *filename, int *istart, int *iend)
{
        ofstream fout (filename) ;
        if(*istart==0)fout << B << " " << volex0 << " " << nxn <<  endl ;
        for(double ib=*istart; ib<=*iend; ib++){
        for(double iq=0; iq<nxn; iq++){
        for(double is=0; is<nxn; is++){
          double mu_b = -muBmx + ib*2*muBmx/(nxn-1.) ;
          double mu_q = -muQmx + iq*2*muQmx/(nxn-1.) ;
          double mu_s = -muSmx + is*2*muSmx/(nxn-1.) ; 
          double T_c = getT_c(mu_b, mu_q, mu_s) ;
          fout << setprecision (7)  << T_c <<  " " ;
          if(is==nxn/2&&iq==nxn/2)cout << " maketctable "<< ib <<" " << setprecision (7) << T_c <<endl;
        } fout  << endl ;
        }}
}

void getmixf(double X, double *fX, double *fpX) 
{
        double Z = X / (1+X/aaa) ;
        double Zp = 1 / ( (1+X/aaa) * (1+X/aaa) );
        *fX  =  exp(-Z-bbb*Z*Z) ;
        *fpX =   (-1-2*bbb*Z) * exp(-Z-bbb*Z*Z) * Zp ;
}


void mix_params(double T, double mu_b, double mu_q, double mu_s, double& T_c , double& dT, double& lam, double& lam_T, double& lam_mub)
{      
        //KW changes to avoid NaN or Inf for large mu

        T_c = tc->pi(mu_b, mu_q, mu_s) ;
        dT = 0.003; //0.00005 ;
        double dlim=3.;        
        double d=mu_b/muc ;                      
        d=min(d,dlim) ;      
        d=max(d,-dlim);                    
        const double del = delta0*exp(-d*d)  ;  
        double del_mub = -2*del *d / muc;       
        if(d>=dlim)del_mub = 0.;                
        double X = (T-T_c)/del ;
        double fX, fpX; 
        getmixf(X,&fX,&fpX);
        if(T>T_c){
                lam  = fX ;
                lam_T = fpX / del  ;
                lam_mub = -fpX * X / del * del_mub ;
                //cout<<X<<" "<<lam<<" "<<lam_T<<" "<<lam_mub<<endl;
        }else{
                lam=1;
                lam_T=0;
                lam_mub=0;
        }
}

void mix_params_old(double T, double mu_b, double mu_q, double mu_s, double& T_c , double& dT, double& lam, double& lam_T, double& lam_mub)
{      
        T_c = tc->pi(mu_b, mu_q, mu_s) ;
        //double T_cc = getT_c(mu_b, mu_q, mu_s) ;
        //cout << " T_c = "<<T_c<< " " << T_cc <<  " mu = " << mu_b<<" "<< mu_q<<" "<< mu_s <<endl;
        dT = 0.00005 ;
        double  alp = 0.5 ; 
        const double delta = delta0*exp(-mu_b*mu_b/muc/muc) ;
        const double delta_mub = -2*delta *mu_b/muc/muc;
        double del = delta * ( 1 + alp * (T-T_c)/T_c ) ;
        double del_T = delta * alp / T_c ;
        double del_mub =  delta_mub * ( 1 + alp * (T-T_c)/T_c ) ;
        if(T>T_c){
                lam  = exp(-(T-T_c)/del);
                lam_T = -lam * ( 1/del - (T-T_c)/del/del*del_T ) ;
                lam_mub = lam*(T-T_c)/del/del*del_mub;
        }else{
                lam=1;
                lam_T=0;
                lam_mub=0;
        }
}

void mix_3f(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{      
        //KW  several changes for better stability
        
        if(T<=0.){
                e = n_b = n_q = n_s = p = 0. ;
                return ;
        }
        double eH, nbH, nqH, nsH, pH, eQ, nbQ, nqQ, nsQ, pQ, lam, lam_T, lam_mub, T_c, dT, tt, ttt[2], aa[5], a[5][2] ;
        int i, ica ;
        mix_params(T, mu_b, mu_q, mu_s, T_c, dT, lam, lam_T, lam_mub);

        if(T<T_c-dT || T>T_c+dT){ica=1; ttt[0]=T;} //exact 
        else {ica=2; ttt[0]=T_c-dT; ttt[1]=T_c+dT;} //interpolation

        for(i=0; i<ica ; i++){
          tt=ttt[i] ;
          HG_exv2(tt, mu_b, mu_q, mu_s, eH, nbH, nqH, nsH, pH) ;
          QGP(tt, mu_b, mu_q, mu_s, eQ, nbQ, nqQ, nsQ, pQ) ;
          mix_params(tt, mu_b, mu_q, mu_s, T_c, dT, lam, lam_T, lam_mub);
          p   = pQ *(1.-lam) + lam*pH ;                                   
          n_b = nbQ*(1.-lam) + lam*nbH  ;                                 
          if(lam_mub!=0.)n_b += lam_mub * (pH-pQ)  ;                      
          n_q = nqQ*(1.-lam) + lam*nqH ;
          n_s = nsQ*(1.-lam) + lam*nsH ;
          e = eQ*(1.-lam)  + lam*eH  ;                                   
          if(tt*lam_T+mu_b*lam_mub!=0.)e+=(tt*lam_T+mu_b*lam_mub)*(pH-pQ);  
          //cout<<"mix_3f: " << tt <<"  "<< mu_b <<"  "<< mu_q <<"  "<< mu_s <<
          //"  ==>  " << eH <<"  "<< nbH <<"  "<< nqH <<"  "<< nsH <<endl;
          if(ica==2){
          a[0][i] = e   ;                           
          a[1][i] = n_b ;                                
          a[2][i] = n_q ; 
          a[3][i] = n_s ; 
          a[4][i] = p   ;                        
          }
        }
        if(ica==2){
          for(int j=0;j<5; j++){
            aa[j] = a[j][0] + (T-T_c+dT)/2./dT * (a[j][1]-a[j][0]) ;
          }
          e    = aa[0] ;                        
          n_b  = aa[1] ;                             
          n_q  = aa[2] ; 
          n_s  = aa[3] ; 
          p    = aa[4] ;                    
        }  

        if(p<0.) p = 0. ;
        if(e<=0.) e = p = n_b = n_q = n_s = 0. ;
}


void mix_simple(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{
        if(T<=0.){
                e = n_b = n_q = n_s = p = 0. ;
                return ;
        }
        double eH, nbH, nqH, nsH, pH, eQ, nbQ, nqQ, nsQ, pQ ;
        const double delta = delta0*exp(-mu_b*mu_b/muc/muc) ;
        HG_exv2(T, mu_b, mu_q, mu_s, eH, nbH, nqH, nsH, pH) ;
        QGP(T, mu_b, mu_q, mu_s, eQ, nbQ, nqQ, nsQ, pQ) ;
        const double lam = 0.5*(1. - (pQ - pH)/sqrt((pQ-pH)*(pQ-pH)+4.*delta)) ;
        p = lam*pH + (1.-lam)*pQ + 2.*delta/sqrt((pQ-pH)*(pQ-pH)+4.*delta) ;
        e = lam*eH + (1.-lam)*eQ - 2.*delta*(1.+mu_b*mu_b/muc/muc)/sqrt((pQ-pH)*(pQ-pH)+4.*delta) ;
        n_b = lam*nbH + (1.-lam)*nbQ - 2.*delta*(mu_b*mu_b/muc/muc)/sqrt((pQ-pH)*(pQ-pH)+4.*delta) ;
        n_q = lam*nqH + (1.-lam)*nqQ ;
        n_s = lam*nsH + (1.-lam)*nsQ ;
        if(p<0.) p = 0. ;
        if(e<=0.) e = p = n_b = n_q = n_s = 0. ;
}




void mix(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p)
{
        if(eosChoice==EOS3F)
        mix_3f(T, mu_b, mu_q, mu_s, e, n_b, n_q, n_s, p) ;
        else{
        mix_simple(T, mu_b, mu_q, mu_s, e, n_b, n_q, n_s, p) ;
        }
}


double fT_HG(double T, double* par)
{
        double &E = par[0] ;
        double &mub = par[1] ;
        double &muq = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        HG(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return e - E ;
}


double fmub_HG(double mub, double* par)
{
        double &NB = par[0] ;
        double &T = par[1] ;
        double &muq = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        HG(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return nb - NB ;
}


double fmuq_HG(double muq, double* par)
{
        double &NQ = par[0] ;
        double &T = par[1] ;
        double &mub = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        HG(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return nq - NQ ;
}


double fmus_HG(double mus, double* par)
{
        double &NS = par[0] ;
        double &T = par[1] ;
        double &mub = par[2] ;
        double &muq = par[3] ;
        double e, nb, nq, ns, p ;
        HG(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return ns - NS ;
}


double fT_QGP(double T, double* par)
{
        double &E = par[0] ;
        double &mub = par[1] ;
        double &muq = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        QGP(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return e - E ;
}


double fmub_QGP(double mub, double* par)
{
        double &NB = par[0] ;
        double &T = par[1] ;
        double &muq = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        QGP(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return nb - NB ;
}


double fmuq_QGP(double muq, double* par)
{
        double &NQ = par[0] ;
        double &T = par[1] ;
        double &mub = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        QGP(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return nq - NQ ;
}


double fmus_QGP(double mus, double* par)
{
        double &NS = par[0] ;
        double &T = par[1] ;
        double &mub = par[2] ;
        double &muq = par[3] ;
        double e, nb, nq, ns, p ;
        QGP(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return ns - NS ;
}


double fT_Mix(double T, double* par)
{
        double &E = par[0] ;
        double &mub = par[1] ;
        double &muq = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        mix(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return e - E ;
}


double fmub_Mix(double mub, double* par)
{
        double &NB = par[0] ;
        double &T = par[1] ;
        double &muq = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        mix(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return nb - NB ;
}


double fmube_Mix(double mub, double* par)
{
        double &E = par[0] ;
        double &T = par[1] ;
        double &muq = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        mix(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return e - E ;
}


double fmuq_Mix(double muq, double* par)
{
        double &NQ = par[0] ;
        double &T = par[1] ;
        double &mub = par[2] ;
        double &mus = par[3] ;
        double e, nb, nq, ns, p ;
        mix(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return nq - NQ ;
}


double fmus_Mix(double mus, double* par)
{
        double &NS = par[0] ;
        double &T = par[1] ;
        double &mub = par[2] ;
        double &muq = par[3] ;
        double e, nb, nq, ns, p ;
        mix(T, mub, muq, mus, e, nb, nq, ns, p) ;
        return ns - NS ;
}


double fT_QGP_uds(double T, double* par)
{
        double &E = par[0] ;
        double &muu = par[1] ;
        double &mud = par[2] ;
        double &mus = par[3] ;
        double e, nu, nd, ns, p ;
        QGP_uds(T, muu, mud, mus, e, nu, nd, ns, p) ;
        return e - E ;
}


double fmuu_QGP_uds(double muu, double* par)
{
        double &NU = par[0] ;
        double &T = par[1] ;
        double &mud = par[2] ;
        double &mus = par[3] ;
        double e, nu, nd, ns, p ;
        QGP_uds(T, muu, mud, mus, e, nu, nd, ns, p) ;
        return nu - NU ;
}

double fmud_QGP_uds(double mud, double* par)
{
        double &ND = par[0] ;
        double &T = par[1] ;
        double &muu = par[2] ;
        double &mus = par[3] ;
        double e, nu, nd, ns, p ;
        QGP_uds(T, muu, mud, mus, e, nu, nd, ns, p) ;
        return nd - ND ;
}


double fmus_QGP_uds(double mus, double* par)
{
        double &NS = par[0] ;
        double &T = par[1] ;
        double &muu = par[2] ;
        double &mud = par[3] ;
        double e, nu, nd, ns, p ;
        QGP_uds(T, muu, mud, mus, e, nu, nd, ns, p) ;
        return ns - NS ;
}


void solveHG(double e, double nb, double nq, double ns, double T, double mub, double muq, double mus)
{
        double mub0, muq0, mus0 ;
        int status ;

        mub = 0. ; muq = 0. ; mus = 0. ;

        for(int i=0; i<10; i++){
        mub0 = mub ; muq0 = muq ; mus0 = mus ;

        double par [4] = {e, mub0, muq0, mus0} ;
        T = findroot(fT_HG, par, 0.04, 0.165, status) ;

        double par2[4] = {nb, T, muq0, mus0} ;
        mub = findroot(fmub_HG, par2, -0.5, 0.5, status) ;

        double par3[4] = {nq, T, mub0, mus0} ;
        muq = findroot(fmuq_HG, par3, -0.5, 0.5, status) ;

        double par4[4] = {ns, T, mub0, muq0} ;
        mus = findroot(fmus_HG, par4, -0.5, 0.5, status) ;

        cout << "root T = " << T << "  mub = " << mub << "  muq = " << muq << "  mus = " << mus << endl ;
        }
}


void solveQGP(double e, double nb, double nq, double ns, double T, double mub, double muq, double mus)
{
        double mub0, muq0, mus0;
        int status ;

        mub = 0. ; muq = 0. ; mus = 0. ;

        if(nb>=nq && nb>=ns){ // baryon potential is principal
        for(int i=0; i<10; i++){
        mub0 = mub ; muq0 = muq ; mus0 = mus ;

        double par [4] = {e, mub0, muq0, mus0} ;
        T = findroot(fT_QGP, par, 0.07, 0.3, status) ;

        double par2[4] = {nb, T, muq, mus} ;
        mub = findroot(fmub_QGP, par2, -0.5, 0.5, status) ;

        double par3[4] = {nq, T, mub, mus} ;
        muq = findroot(fmuq_QGP, par3, -0.5, 0.5, status) ;

        double par4[4] = {ns, T, mub, muq} ;
        mus = findroot(fmus_QGP, par4, -0.5, 0.5, status) ;

        cout << "root T = " << T << "  mub = " << mub << "  muq = " << muq << "  mus = " << mus << endl ;
        }
        return ;
        }

        if(nq>=nb && nq>=ns){ // charge chem. potential is principal
        for(int i=0; i<10; i++){
        mub0 = mub ; muq0 = muq ; mus0 = mus ;

        double par [4] = {e, mub0, muq0, mus0} ;
        T = findroot(fT_QGP, par, 0.07, 0.3, status) ;

        double par3[4] = {nq, T, mub, mus} ;
        muq = findroot(fmuq_QGP, par3, -0.5, 0.5, status) ;

        double par2[4] = {nb, T, muq, mus} ;
        mub = findroot(fmub_QGP, par2, -0.5, 0.5, status) ;

        double par4[4] = {ns, T, mub, muq} ;
        mus = findroot(fmus_QGP, par4, -0.5, 0.5, status) ;

        cout << "root T = " << T << "  mub = " << mub << "  muq = " << muq << "  mus = " << mus << endl ;
        }
        return ;
        }

        if(ns>=nb && ns>=nq){ // strange chem. potential is principal
        for(int i=0; i<10; i++){
        mub0 = mub ; muq0 = muq ; mus0 = mus ;

        double par [4] = {e, mub0, muq0, mus0} ;
        T = findroot(fT_QGP, par, 0.07, 0.3, status) ;

        double par4[4] = {ns, T, mub, muq} ;
        mus = findroot(fmus_QGP, par4, -0.5, 0.5, status) ;

        double par3[4] = {nq, T, mub, mus} ;
        muq = findroot(fmuq_QGP, par3, -0.5, 0.5, status) ;

        double par2[4] = {nb, T, muq, mus} ;
        mub = findroot(fmub_QGP, par2, -0.5, 0.5, status) ;

        cout << "root T = " << T << "  mub = " << mub << "  muq = " << muq << "  mus = " << mus << endl ;
        }
         return ;
        }
}


void solveHG_lowe(double e, double nb, double nq, double ns, double &T, double &mub, double &muq, double &mus, double &p)
{
        double par [4] = {e, 0., 0., 0.} ;
        int status ;
        T = findroot(fT_HG, par, 0., 0.3, status) ;
        mub = muq = mus = 0. ;
        HG(T, 0., 0., 0., e, nb, nq, ns, p) ;
}

void setaccur(int i, double T, double mub, double muq, double mus)
{
        if(i>0 && fabs(T-tc->pi(mub, muq, mus))<0.001){Raccuracy=1e-7 ; Taccuracy = 0.000005;}
        else {Raccuracy=1e-4 ; Taccuracy = 0.0001;}
}

//---------------------------------------------------------------------------------------------------
   void selectroot(int i, double e, double nb, double nq, double ns, 
   double T , double mub , double muq , double mus , int j, double &err1)
//---------------------------------------------------------------------------------------------------
{
        //Kw Error computing (basis for root selection)

        int idebug = 0 ;
       
        double  exx, nbxx, nqxx, nsxx, p;
        mix(T, mub, muq, mus, exx, nbxx, nqxx, nsxx, p) ;  
        err1 = fabs(exx-e)/max(0.01,e) + fabs(nbxx-nb)/max(0.01,nb) 
                  + fabs(nqxx-nq)/max(0.01,nq) + fabs(nsxx-ns)/max(0.01,ns) ;

        if(idebug)cout <<"    ====> "<<i<<"   " << e <<" "<< nb <<" "
          << nq <<" "<< ns <<"    " << T<<" "<< mub <<" "<< muq <<" "<< mus <<"   ROOT"<<j+1<<"   "<<err1<< endl;  

        //or:
        //cout <<"         => "<<T0<<" "<<T<<" "<<mub<<" "<<muq<<" "<<mus<<"  =>  ";
        //cout            << e <<" "<<exx<<"  "<< nb <<" "<<nbxx<<"  "<< nq <<" "<<nqxx<<"  "<< ns <<" "<<nsxx <<endl  ;   

}
//---------------------------------------------------------------------------------------------------
  void findroot2(int j, double (*f)(double, double*), double* par, double xmin, double xmax, int &nroo,double* rts)
//---------------------------------------------------------------------------------------------------
{
	if(f(xmin, par)*f(xmax, par)>0.){nroo=1;
		if(fabs(f(xmin, par))<fabs(f(xmax, par)))  rts[0] = xmin ;
		else rts[0] =  xmax ;
                return;
	}

        //KW Root bracketing

        int idebug = 0 ;

        // j: 1 = T | mub,muq,mus   2 = mub | T,muq,mus   3 = muq | T,mub,mus   4 = mus | T,mub,muq  

        double delx=0.2 ,x1=xmax,x2=x1-delx,ab[10]   ;
        x2=min(2.,x2); 
        int i = 0;
        
        //cout <<"--++-- "<<xmax*10<<"  "<< f(xmax*10, par) <<endl;
        while(x1>xmin){
           if(idebug)cout <<"----"<<j<<"---- "<<setw(12)<<x1<<"  "<< f(x1, par)+par[0] <<"      "<< par[0]<<"      " << par[1]<<"  "<< par[2]<<"  "<< par[3]<<" "  <<endl;
           if(f(x1, par)*f(x2, par)<0.){
               if(i==10)cout<<"ERROR: More than 5 roots!!!"<<endl;
               else{ 
               i++;ab[i-1]=x2;
               i++;ab[i-1]=x1;
               if(idebug)cout<<"----"<<j<<"---- bracket: "<<x2<<" "<<x1<<"      "<<endl;
               }
           }
           x1=x2;x2=x1-delx;
           if(x2<-2.&&xmin<-2)x2=xmin; 
           if(x2<xmin&&x1>xmin)x2=xmin;
        }
        if(idebug)cout <<"----"<<j<<"---- "<<setw(12)<<x1 <<"  "<<f(x1, par)+par[0] <<endl;
        //cout <<"--++-- "<<xmin*10<<"  "<< f(xmin*10, par) <<endl;

        nroo=i/2;
        if(nroo==0)return; 

        //KW Root finding (up to 3)

        for(i=0;i<nroo;i++){
           xmin=ab[2*i];
           xmax=ab[2*i+1];
	   do{
           	if(f((xmin+xmax)/2., par)*f(xmin, par)>0.) xmin = (xmin+xmax)/2. ;
                else xmax = (xmin+xmax)/2. ;   
                //cout<<xmin<<"  "<<f(xmin, par)+par[0]<<"    "<<xmax<<"  "<<f(xmax, par)+par[0]<<endl;
           }while(fabs(xmax-xmin)>Raccuracy) ;
                //cout<< (xmin+xmax)/2. << " " << fabs(xmax-xmin) <<endl;
           rts[i]=(xmin+xmax)/2. ;
        }
}
//---------------------------------------------------------------------------------------------------------------------
  void SolveMix2(double e, double nb, double nq, double ns, double &T, double &mub, double &muq, double &mus, double &p)
//---------------------------------------------------------------------------------------------------------------------
{
        //KW Allowing non-monotonic functions

        //int cas=0;
        int j, jx;
        double  exx, nbxx, nqxx, nsxx, rts [5], err1, err1x ;                 
        if(e<=0.){
                T = mub = muq = mus = p = 0. ;
                return ;
        }
        double mub0, muq0, mus0, T0;

        mub = 0. ; muq = 0. ; mus = 0. ; T = 0. ;
        int i = 0 ,nroo=0;

        if(fabs(nb)>=fabs(nq) && fabs(nb)>=fabs(ns)){ // baryon potential is principal
        do{
        T0 = T ;
        mub0 = mub ; muq0 = muq ; mus0 = mus ;
 
        double par [4] = {e, mub0, muq0, mus0} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(1, fT_Mix, par, 0., temmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){T=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        T=rts[jx];
        }
        double par2[4] = {nb, T, muq, mus} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(2, fmub_Mix, par2, -muBmax, muBmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){mub=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        mub=rts[jx];
        }
        //if(nroo>1){cout<<"CHOSEN "<<jx<<" "<<mub <<endl;}
        double par3[4] = {nq, T, mub, mus} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(3, fmuq_Mix, par3, -muQmax, muQmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){muq=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        muq=rts[jx];
        }
        double par4[4] = {ns, T, mub, muq} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(4, fmus_Mix, par4, -muSmax, muSmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){mus=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        mus=rts[jx];
        }
        setaccur(i, T, mub, muq, mus) ;
        mix(T, mub, muq, mus, exx, nbxx, nqxx, nsxx, p) ;
        if(exx>0.) p = p * e/exx ;  //new 3118
        i++ ;
        }while(i<maxIter /* && fabs(mub)<muBmax && fabs(muq)<muQmax &&  fabs(mus)<muSmax*/ 
                 && ( fabs(T-T0)/T>Taccuracy || fabs(exx-e)/e>0.01 || fabs(nbxx-nb)>0.01 
                       || fabs(nqxx-nq)>0.01 || fabs(nsxx-ns)>0.01 ) ) ;
        //cas=1;
        goto End ;
        }

        if(fabs(nq)>=fabs(nb) && fabs(nq)>=fabs(ns)){ // charge chem. potential is principal
        do{
        T0 = T ; mub0 = mub ; muq0 = muq ; mus0 = mus ;

        double par [4] = {e, mub0, muq0, mus0} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(1, fT_Mix, par, 0., temmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){T=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        T=rts[jx];
        }
        double par3[4] = {nq, T, mub, mus} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(3, fmuq_Mix, par3, -muQmax, muQmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){muq=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        muq=rts[jx];
        }
        double par2[4] = {nb, T, muq, mus} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(2, fmub_Mix, par2, -muBmax, muBmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){mub=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        mub=rts[jx];
        }
        double par4[4] = {ns, T, mub, muq} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(4, fmus_Mix, par4, -muSmax, muSmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){mus=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        mus=rts[jx];
        }
        setaccur(i, T, mub, muq, mus) ;
        mix(T, mub, muq, mus, exx, nbxx, nqxx, nsxx, p) ;
        if(exx>0.) p = p * e/exx ;  //new 3118
        i++ ;
        }while(i<maxIter && fabs(mub)<muBmax && fabs(muq)<muQmax &&  fabs(mus)<muSmax 
               && ( fabs(T-T0)/T>Taccuracy || fabs(exx-e)/e>0.01 || fabs(nbxx-nb)>0.01 
                    || fabs(nqxx-nq)>0.01 || fabs(nsxx-ns)>0.01 ) ) ;
        //cas=2;
        goto End ;
        }

        if(fabs(ns)>=fabs(nb) && fabs(ns)>=fabs(nq)){ // strange chem. potential is principal
        do{
        T0 = T ; mub0 = mub ; muq0 = muq ; mus0 = mus ;

        double par [4] = {e, mub0, muq0, mus0} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(1, fT_Mix, par, 0., temmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){T=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        T=rts[jx];
        }
        double par4[4] = {ns, T, mub, muq} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(4, fmus_Mix, par4, -muSmax, muSmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){mus=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        mus=rts[jx];
        }
        double par3[4] = {nq, T, mub, mus} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(3, fmuq_Mix, par3, -muQmax, muQmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){muq=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        muq=rts[jx];
        }
        double par2[4] = {nb, T, muq, mus} ;
        setaccur(i, T, mub, muq, mus) ;
        findroot2(2, fmub_Mix, par2, -muBmax, muBmax, nroo, rts) ;
        if(nroo==0){goto End;} 
        for(j=0,jx=0,err1x=1e100;j<nroo;j++){mub=rts[j] ;
        selectroot(i, e, nb, nq, ns, T, mub ,muq, mus , j, err1); 
        if(err1<err1x){err1x=err1;jx=j;}
        mub=rts[jx];
        }
        setaccur(i, T, mub, muq, mus) ;
        mix(T, mub, muq, mus, exx, nbxx, nqxx, nsxx, p) ;
        if(exx>0.) p = p * e/exx ;  //new 3118
        i++ ;
        }while(i<maxIter && fabs(mub)<muBmax && fabs(muq)<muQmax &&  fabs(mus)<muSmax 
               && ( fabs(T-T0)/T>Taccuracy || fabs(exx-e)/e>0.01 || fabs(nbxx-nb)>0.01 
                    || fabs(nqxx-nq)>0.01 || fabs(nsxx-ns)>0.01 ) ) ;
        //cas=3;
        goto End ;
        }
        
        End:
        if(nroo==0)
        {
        T=999.;mub=999.;muq=999.;mus=999.;p=999.;
        }
        else{//only printout for testing
        //double ex, nbx, nqx, nsx, px;
        //mix(T, mub, muq, mus, ex, nbx, nqx, nsx, px) ;
        //cout <<"OK"<< cas<<" =======> "<< T <<" "<< mub  <<"   "<< muq  <<" "<< mus  
        //     <<"  =========> "<< e <<" "<< ex  <<"   "<< nb  <<" "<< nbx  
        //     <<"   "<< nq <<" "<< nqx    <<"   "<< ns <<" "<< nsx<<endl;  
        } 
}//------------------------------------------------------------------------------------------------------------------------

double findroot(double (*f)(double, double*), double* par, double xmin, double xmax, int &status)
{
        if(f(xmin, par)*f(xmax, par)>0.){
                status = S_NOROOT ;
                if(fabs(f(xmin, par))<fabs(f(xmax, par))) return xmin ;
                else return xmax ;
        }
	do{
           	if(f((xmin+xmax)/2., par)*f(xmin, par)>0.) xmin = (xmin+xmax)/2. ;
                else xmax = (xmin+xmax)/2. ;
        }while(fabs(xmax-xmin)>Raccuracy) ;
                status = S_OK ;
                //cout<< (xmin+xmax)/2. << " " << fabs(xmax-xmin) <<endl;
                return (xmin+xmax)/2. ;
}

double deltaMu(double e, double nb, double T, int& status, int& status2)
{
        double par [4] = {e, T, 0., 0.} ;
        double mub_e = findroot(fmube_Mix, par, 0., 1., status) ;
        double parn [4] = {nb, T, 0., 0.} ;
        double mub_n = findroot(fmub_Mix, parn, 0., 1., status2) ;
        return mub_e - mub_n ;
}


// solves {e,nb} --> {T, mub} relations
void solve2(double e, double nb, double &T, double &mu, double &p)
{
        int status, status2 ;
        double nq, ns ;
        double nsign = nb!=0. ? nb/fabs(nb) : 0. ;
        nb = fabs(nb) ;
//        cout << "-----\n" ;
//        for(double tt=0.149; tt<0.151; tt+=0.0001)
//                cout << setw(14) << tt << setw(14) << deltaMu(e, nb, tt, status, status2) << endl ;
//        cout << "-----\n" ;
        double parTmax [4] = {e, 0., 0., 0.} ;
        double Tmax = findroot(fT_Mix, parTmax, 0., 0.5, status) ;
//        cout << "Tmax = " << Tmax << endl ;
        double parTmin [4] = {e, muBmax, 0., 0.} ;
        double Tmin = findroot(fT_Mix, parTmin, 0., 0.5, status) ;
//        cout << "Tmin = " << Tmin << endl ;
//---- extreme cases
        double dmu1 = deltaMu(e, nb, Tmin, status, status2) ;
        double dmu2 = deltaMu(e, nb, Tmax, status, status2) ;
//        cout << "dmu1= " << dmu1 << " dmu2= " << dmu2 << endl ;
        if(dmu1>0. && dmu2>0.){
                T = Tmax ; mu = 0. ; 
                mix(T, mu, 0., 0., e, nb, nq, ns, p) ;
                return ;
        }
        if(dmu1<0. && dmu2<0.){
                T = Tmin ; mu = muBmax ;
                mix(T, mu, 0., 0., e, nb, nq, ns, p) ;
                return ;
        }
//--------
        double dmu ;
        do{
                T = (Tmin+Tmax)/2. ;
                dmu = deltaMu(e, nb, T, status, status2) ;
                if(dmu>0.) Tmin = T ;
                else Tmax = T ;
        }while((Tmax-Tmin)/Tmax>0.0001) ;
        dmu1 = deltaMu(e, nb, Tmin, status, status2) ;
        dmu2 = deltaMu(e, nb, Tmax, status, status2) ;
        if(fabs(dmu1)>fabs(dmu2)) T = Tmax ; else  T = Tmin ;
        double parn [4] = {nb, T, 0., 0.} ;
        mu = findroot(fmub_Mix, parn, 0., muBmax, status2) ;
        mu *= nsign ;
        mix(T, mu, 0., 0., e, nb, nq, ns, p) ;

//        exit(1) ;
//        for(double T=0.; T<Tmax; T+=0.002){
//        double par [4] = {e, T, 0., 0.} ;
//        double mub_e = findroot(fmube_Mix, par, 0., 1., status) ;
//        double parn [4] = {nb, T, 0., 0.} ;
//        double mub_n = findroot(fmub_Mix, parn, 0., 1., status2) ;
//        if(status==S_OK && status2==S_OK) cout << " mu_e = " << setw(14) << mub_e << " mu_n = " << setw(14) << mub_n << 
//                " T = " << setw(5) << T << " delta = " << setw(14) << mub_e-mub_n << endl ;
//        }
}


void solveQGP_uds(double e, double nu, double nd, double ns, double &T, double &muu, double &mud, double &mus, double &p)
{
        if(e<=0.){
                T = muu = mud = mus = p = 0. ;
                return ;
        }
        if(e<0.5){
                solveHG_lowe(e, 0., 0., 0., T, muu, mud, mus, p) ;
                return ;
        }
        double muu0, mud0, mus0;
        int status ;

        muu = 0. ; mud = 0. ; mus = 0. ;

        if(fabs(nu)>=fabs(nd) && fabs(nu)>=fabs(ns)){ // baryon potential is principal
        for(int i=0; i<4; i++){
        muu0 = muu ; mud0 = mud ; mus0 = mus ;

        double par [4] = {e, muu0, mud0, mus0} ;
        T = findroot(fT_QGP_uds, par, 0.04, 0.5, status) ;

        double par2[4] = {nu, T, mud, mus} ;
        muu = findroot(fmuu_QGP_uds, par2, -1.0, 1.0, status) ;

        double par3[4] = {nd, T, muu, mus} ;
        mud = findroot(fmud_QGP_uds, par3, -1.0, 1.0, status) ;

        double par4[4] = {ns, T, muu, mud} ;
        mus = findroot(fmus_QGP_uds, par4, -1.0, 1.0, status) ;

//        cout << "root T = " << T << "  muu = " << muu << "  mud = " << mud << "  mus = " << mus << endl ;
        }
        QGP_uds(T, muu, mud, mus, e, nu, nd, ns, p) ;
        return ;
        }

        if(fabs(nd)>=fabs(nu) && fabs(nd)>=fabs(ns)){ // charge chem. potential is principal
        for(int i=0; i<4; i++){
        muu0 = muu ; mud0 = mud ; mus0 = mus ;

        double par [4] = {e, muu0, mud0, mus0} ;
        T = findroot(fT_QGP_uds, par, 0.04, 0.5, status) ;

        double par3[4] = {nd, T, muu, mus} ;
        mud = findroot(fmud_QGP_uds, par3, -1.0, 1.0, status) ;

        double par2[4] = {nu, T, mud, mus} ;
        muu = findroot(fmuu_QGP_uds, par2, -1.0, 1.0, status) ;

        double par4[4] = {ns, T, muu, mud} ;
        mus = findroot(fmus_QGP_uds, par4, -1.0, 1.0, status) ;

//        cout << "root T = " << T << "  muu = " << muu << "  mud = " << mud << "  mus = " << mus << endl ;
        }
        QGP_uds(T, muu, mud, mus, e, nu, nd, ns, p) ;
        return ;
        }

        if(fabs(ns)>=fabs(nu) && fabs(ns)>=fabs(nd)){ // strange chem. potential is principal
        for(int i=0; i<4; i++){
        muu0 = muu ; mud0 = mud ; mus0 = mus ;
        
        double par [4] = {e, muu0, mud0, mus0} ;
        T = findroot(fT_QGP_uds, par, 0.04, 0.5, status) ;

        double par4[4] = {ns, T, muu, mud} ;
        mus = findroot(fmus_QGP_uds, par4, -1.0, 1.0, status) ;

        double par3[4] = {nd, T, muu, mus} ;
        mud = findroot(fmud_QGP_uds, par3, -1.0, 1.0, status) ;

        double par2[4] = {nu, T, mud, mus} ;
        muu = findroot(fmuu_QGP_uds, par2, -1.0, 1.0, status) ;

//        cout << "root T = " << T << "  muu = " << muu << "  mud = " << mud << "  mus = " << mus << endl ;
        }
        QGP_uds(T, muu, mud, mus, e, nu, nd, ns, p) ;
        return ;
        }
}


// generates a part of 1f table
void makeInvTable2(const char *filename, int istart, int iend)
{
        ofstream fout (filename) ;
        if(istart==0)
        fout << setw(15) << emax << setw(15) << e0 
                <<  setw(15) << nxe << setw(15) << nxn << endl ;
        double T, mub, p ;
        for(double ixe=istart; ixe<iend; ixe++)
        for(double ixnb=0; ixnb<nxn; ixnb++){
                double xe = ixe*del_xe ;
                double e = e0*(exp(xe)-1.) ;
                double nmax = anmax * pow(e,bnmax); 
                double n0 = an0*pow(0.5*e, bnmax) ;
                double rat = 0 ;
                if(n0>0.0) rat = nmax/n0 ; 
                double xn_max = log(rat+1) ;
                double del_xn = 2.*xn_max/(nxn-1) ;
                double xnb = -xn_max + ixnb*del_xn ;
                double nb = xnb/fabs(xnb+1e-20)*n0*(exp(fabs(xnb))-1.) ;
                solve2(e, nb, T, mub, p) ;
                fout << p << " " << T << " " << mub << endl ;
                cout << "ixnb = " << ixnb << endl ;
        }
        fout.close() ;
}


//---------------------------------------------------------------------------------------------------
  void makeinvtable_(char *fnamemtr, char *filename, int *istart, int *iend, int *jstart, int *jend)
//---------------------------------------------------------------------------------------------------
{
        if(strlen(fnamemtr)>0) {
                m2file.open(fnamemtr,ios_base::app) ;
                cout.rdbuf(m2file.rdbuf()) ;
        }
        cout <<"makeinvtable -- mtr:"<< fnamemtr<<"  out:"<< filename<< "  istart:"<< *istart<< "  iend:"<< *iend<< "  jstart:"<< *jstart<< "  jend:"<< *jend<<endl;
        ofstream fout (filename,ios::app) ;
        //-----generating e / n grids
        double egrid [nxe], ngrid[nxe][nxn];
        for(int ixe=0; ixe<nxe; ixe++){
                double xe = ixe * del_xe ;
                egrid[ixe] = e0 * (exp(xe)-1.);
                double nmax = anmax * pow(egrid[ixe],bnmax);  // <----- nmax as a function of e !!
                double n0 = an0*pow(0.5*egrid[ixe], bnmax) ;
                double rat = 0 ;
                if(n0>0.0) rat = nmax/n0 ; 
                //cout<<"makeinvtable: "<<*istart<<" "<<*iend<<"  "<<ixe<<" "<<egrid[ixe]<<" "<<rat<<endl;
                for(int ixn=0; ixn<nxn; ixn++){
                        double xn_max = log(rat+1) ;
                        double del_xn = 2.*xn_max/(nxn-1) ;
                        double xn = -xn_max + ixn*del_xn ;
                        ngrid[ixe][ixn] =  xn/fabs(xn+1e-20)*n0*(exp(fabs(xn))-1.) ;
                }
        }
        //----writing grid to file if istart=jstart=0 
        if(*istart==0 && *jstart==0){
        fout << "  " << B << "  " << volex0 << "  " << delta0 << "  " << aaa <<  "  " << bbb <<  endl ;
        fout << "  " << emax << "  " << e0   <<  "  " << nxe << "  " << nxn << endl ;
        for(int ixe=0;  ixe<nxe;  ixe++)
                fout << egrid[ixe] << " "  ;
        fout << endl ;
        for(int ixe=0; ixe<nxe; ixe++)
        for(int ixn=0; ixn<nxn; ixn++)
                fout << ngrid[ixe][ixn] << " "  ;     
        fout << endl ;
        }
        //------------------------------------------    
        double T, mub, muq, mus, p ;
        for(int ixe =*istart; ixe <=*iend; ixe++)
        for(int ixnb=*jstart; ixnb<=*jend; ixnb++)
        for(int ixnq=0; ixnq<nxn; ixnq++)
        for(int ixns=0; ixns<nxn; ixns++){
                if(ixns==0) cout << "   e(" <<ixe<< ") = " << egrid[ixe] << "   nb("
                        << ixnb << ") = " << ngrid[ixe][ixnb] << "   nq("<< ixnq << ") = " << ngrid[ixe][ixnq] << endl ;
                SolveMix2(egrid[ixe], ngrid[ixe][ixnb], ngrid[ixe][ixnq], ngrid[ixe][ixns], T, mub, muq, mus, p) ;
                fout  << p << " " << T << " " << mub << " " << muq << " " << mus << endl ;
                //cout << temmax<< " In  e nb nq ns :     " << e << " " << nb << " " << nq << " " << ns  << endl ;
                //cout << temmax<< " Out p T mub muq mus :" << p << " " << T << " " << mub << " " << muq << " " << mus << endl ;
                //fout     << ixe << " " << ixnb << " " << ixnq << " " << ixns << "         " << e << " " << nb << " " << nq << " " << ns << "         " << p << " " << T << " " << mub << " " << muq << " " << mus << endl ;
        }
        //---------------------------------------- 
        fout.close() ;
}//-----------------------------------------------------------------------------


void makeinvtable1f_(char *fnamemtr, char *filename, int *istart, int *iend, int *jstart, int *jend)
{
 if(strlen(fnamemtr)>0) {
 m2file.open(fnamemtr,ios_base::app) ;
 cout.rdbuf(m2file.rdbuf()) ;
 }
 cout <<"makeinvtable -- mtr:"<< fnamemtr<<"  out:"<< filename<< "  istart:"<< *istart<< "  iend:"<< *iend
 << "  jstart:"<< *jstart<< "  jend:"<< *jend<<endl;
 ofstream fout (filename,ios::app) ;
 //-----generating e / n grids
 double egrid [nxe], ngrid[nxe][nxn];
        for(int ixe=0; ixe<nxe; ixe++){
                double xe = ixe * del_xe ;
         egrid[ixe] = e0 * (exp(xe)-1.);
  double nmax = anmax * pow(egrid[ixe],bnmax); 
  double n0 = an0*pow(0.5*egrid[ixe], bnmax) ;
  double rat = 0 ;
  if(n0>0.0) rat = nmax/n0 ; 
  for(int ixn=0; ixn<nxn; ixn++){
   double xn_max = log(rat+1) ;
   double del_xn = 2.*xn_max/(nxn-1) ;
                 double xn = -xn_max + ixn*del_xn ;
                 ngrid[ixe][ixn] =  xn/fabs(xn+1e-20)*n0*(exp(fabs(xn))-1.) ;
  }
        }
 //----writing grid to file if istart=jstart=0 
        if(*istart==0 && *jstart==0){
        fout << "  " << B << "  " << volex0 << "  " << delta0 << "  " << aaa <<  "  " << bbb <<  endl ;
        fout << "  " << emax << "  " << e0 
                <<  "  " << nxe << "  " << nxn << endl ;
        for(int ixe=0;  ixe<nxe;  ixe++) { 
                fout << egrid[ixe] << " "  ;
        }
        fout << endl ;
 for(int ixe=0; ixe<nxe; ixe++)
        for(int ixn=0; ixn<nxn; ixn++)
                fout << ngrid[ixe][ixn] << " "  ;     
 fout << endl ;
 }
 //------------------------------------------    
        double nq=0., ns=0., T, mub, muq, mus, p ;
        for(int ixe =*istart; ixe <=*iend; ixe++){
  cout << ixe << " / " << *iend << endl;
        for(int ixnb=*jstart; ixnb<=*jend; ixnb++){
  eosohlle_(&egrid[ixe], &ngrid[ixe][ixnb], &nq, &ns, &T, &mub, &muq, &mus, &p) ;
  fout  << p << " " << T << " " << mub << endl ;
         //cout << temmax<< " In  e nb nq ns :     " << e << " " << nb << " " << nq << " " << ns  << endl ;
         //cout << temmax<< " Out p T mub muq mus :" << p << " " << T << " " << mub << " " << muq << " " << mus << endl ;
  //fout     << ixe << " " << ixnb << " " << ixnq << " " << ixns << "         " << e << " " << nb << " " << nq << " " << ns << "         " << p << " " << T << " " << mub << " " << muq << " " << mus << endl ;
        }}
 //---------------------------------------- 
        fout.close() ;
}


//###################################################################################
//###################################################################################
//###################################################################################
//###################################################################################
//###################################################################################
//###################################################################################
//###################################################################################
//###################################################################################


int igetnparticles_(void)  // all
{
        //return database->GetNParticles() ;
        return database->GetNParticles(true) ;
}

int igetkparticles_(void)   // good + charm
{
        //return database->GetNParticles() ;
        return database->GetKParticles() ;
}

void convertToFString(char *s)
{
        int ilast=0 ;
        while(s[ilast]!=0x0) ilast++ ;
        for(int i=ilast; i<cStringLength; i++) s[i] = ' ' ;
//        s[cStringLength-1] = 0x0 ; // to make it a valid C string, not needed here
}


void igetparticle_(char* name, double* mass, int* charge, int* id_PDG, 
        int* degeneracy, int* stats ,int* baryon_number, int* strangeness, int* icode, int* iflag)
// icode : see description below
// iflag return value:
//         1: particle belongs to EPOS1.99 list,
//         2: particle is in a wider list charm,bottom,top=0 (including EPOS1.99 list)
//         3: otherwise
{
 static int iparticleindex = 0;  
 if(*icode==0)iparticleindex=0;
 const int epos199list [] = {111, 211, 221, 311, 321, 331, 113, 213, 223, 323, 313, 333,
   2212, 2112, 3122, 3222, 3212, 3112, 3322, 3312, 2224, 2214, 2114, 1114, 3224, 3214, 3114,
   3324, 3314, 3334};   // hadrons in EPOS 1.99
 ParticlePDG2* particle = database->GetPDGParticleByIndex(iparticleindex);
 if(!particle) { cout << "Invalid particle index!\n" ; exit(1) ; }
 strcpy(name, particle->GetName()) ;
 *mass = particle->GetMass() ;
 *charge = particle->GetElectricCharge() ;
 *id_PDG = particle->GetPDG() ;
 *degeneracy = int(2.*particle->GetSpin() + 1.) ;
 *stats = int(2.*particle->GetSpin()) & 1 ? 1 : -1 ;
 *baryon_number = particle->GetBaryonNumber() ;
 *strangeness = particle->GetStrangeness() ;
 convertToFString(name) ;
 //KW:
 //icode :  good && !last : 1  ; good && last : -1    <---used in EPOS
 //         bad && !last : 9999  ; bad && last : -9999
 //definition of good/bad status: as before, but charm is allowed
 *icode = 1 ;
 double width = particle->GetWidth() ;
 Char_t staPDG = particle->GetPDGstatus() ;
 if(iparticleindex>0){
   ParticlePDG2* particle1 = database->GetPDGParticleByIndex(iparticleindex-1);
   if(abs(*id_PDG)==abs(particle1->GetPDG()))staPDG=particle1->GetPDGstatus();} //not defined for antibaryons
 int charm = particle->GetCharm() ;
 int bcharge = particle->GetBeauty() ;
 int tcharge = particle->GetTruth() ;
 Double_t fMinimumMass = database->GetMinimumMass() ;
 Double_t fMaximumMass = database->GetMaximumMass() ;
 Double_t fMinimumWidth = database->GetMinimumWidth() ;
 Double_t fMaximumWidth = database->GetMaximumWidth() ;
 bool inEpos199 = false; // flag if the particle is in EPOS 1.99 list
 for(unsigned int i=0; i<sizeof(epos199list)/sizeof(int); i++)
  if(abs(*id_PDG)==epos199list[i]) inEpos199 = true;
 *iflag = 3;
 if(charm==0 && bcharge==0 && tcharge==0) *iflag = 2;
 if(inEpos199) *iflag = 1;
 if(!(fMinimumMass<=*mass && *mass<=fMaximumMass) || !(fMinimumWidth<=width && width<=fMaximumWidth) || abs(*id_PDG)<100
 || staPDG=='S' || staPDG=='F' || bcharge!=0  || tcharge!=0)  *icode = 9999 ;
 iparticleindex++ ;
 if(iparticleindex==database->GetNParticles(true)){
  *icode = *icode*(-1);
  iparticleindex=0;
 }
 //cout << setw(10) << iparticleindex << setw(10) << *id_PDG  << setw(10) << *icode    << setw(10) << *mass << setw(10) << width << setw(10) << staPDG<< setw(10) << bcharge  << endl;
}


