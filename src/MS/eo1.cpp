//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

//#include <TError.h>
//#include <TApplication.h>
//#include <TGraph.h>
//#include <TCanvas.h>
//#include <TMath.h>
//#include <TGraph.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
//#include <TF1.h>

#include "eo1.h"

using namespace std ;


EoS1f::EoS1f(char *filename)
{
	ifstream fin (filename) ;

	fin >> emax >> e0 >> nmax >> n0 >> ne >> nn ;

        egrid = new double [ne] ;
        ngrid = new double* [ne] ;
        for(int i=0; i<ne; i++)
                ngrid[i] = new double [nn] ;        

	for(int ixe=0;  ixe<ne;  ixe++) fin >> egrid[ixe]  ; 
        
        for(int ixe=0;  ixe<ne;  ixe++)
        for(int ixnb=0; ixnb<nn; ixnb++)
                fin >> ngrid[ixe][ixnb] ;
               
	T = new double [ne*nn] ;
	pre = new double [ne*nn] ;
	mub = new double [ne*nn] ;

	for(int ie=0; ie<ne; ie++)
	for(int inb=0; inb<nn; inb++){
		fin >> pre[index(ie,inb)] >> T[index(ie,inb)] 
		>> mub[index(ie,inb)]  ;
                if( pre[index(ie,inb)] < 1e-18 ) pre[index(ie,inb)] = 0. ;
	}
	fin.close() ;
}

EoS1f::~EoS1f()
{
	delete [] T ;
	delete [] pre ;
	delete [] mub ;
}

void EoS1f::eosranges(double &_emax, double &_e0, double &_nmax, double &_n0, int &_ne, int &_nn)
{
        _emax=emax;
        _e0=e0;
        _nmax=nmax;
        _n0=n0;
        _ne=ne;
        _nn=nn;
}

void EoS1f::getue(double e, int &ixe, double &ue, int &iout){
        int k, k1, k2;
        k1 = 0 ;
        k2 = ne-1 ;
        for(;k2-k1>1;){k=(k1+k2)/2 ; if(e<egrid[k])k2=k ; else k1=k ;}
        ixe=k1;
	if(ixe>ne-2) ixe = ne - 2 ;
	double dxe = egrid[k2] - egrid[k1] ;
	double xem = e - egrid[k1] ;
        ue = xem/dxe ;
} 

void EoS1f::getun(int ie, double n, int &ixn_low, int &ixn_up, double &un_low, double &un_up, int &iout){
        int k, k1, k2;
        // --- lower side, so for ngird[ie][...]
        k1 = 0 ;
        k2 = nn-1 ;
        for(;k2-k1>1;){k=(k1+k2)/2 ; if(n<ngrid[ie][k])k2=k ; else k1=k ;}
        ixn_low=k1;
        if(ixn_low>nn-2) ixn_low = nn - 2 ;
        if(ixn_low<0) ixn_low = 0 ;
        double dxn = ngrid[ie][k2] - ngrid[ie][k1] ;
        double xnm=n - ngrid[ie][k1] ;
        un_low = xnm/dxn ;
        // --- upper side so for ngird[ie+1][...]
        k1 = 0 ;
        k2 = nn-1 ;
        for(;k2-k1>1;){k=(k1+k2)/2 ; if(n<ngrid[ie+1][k])k2=k ; else k1=k ;}
        ixn_up=k1;
        if(ixn_up>nn-2) ixn_up = nn - 2 ;
        if(ixn_up<0) ixn_up = 0 ;
        dxn = ngrid[ie+1][k2] - ngrid[ie+1][k1] ;
        xnm=n - ngrid[ie+1][k1] ;
        un_up = xnm/dxn ;
}

void EoS1f::eos(double e, double nb, double nq, double ns,
		double &_T, double &_mub, double &_muq, double &_mus, double &_p)
{
	if(e<0.) { _T = _mub = _muq = _mus = _p = 0. ; }
        int ixe , ixnb_low, ixnb_up, iout;
        double ue, ub_low, ub_up ;
          
        getue(e,  ixe , ue, iout);
        getun(ixe, nb, ixnb_low, ixnb_up, ub_low, ub_up, iout);
               
	double  we [2] = {1.-ue, ue} ;
        const double wnb [2][2] = {{1.-ub_low, ub_low},{1.-ub_up, ub_up}} ;
        const int ixnb[2] = {ixnb_low, ixnb_up} ;

	_T = _mub = _muq = _mus = _p = 0. ;
	for(int je=0; je<2; je++)
	for(int jnb=0; jnb<2; jnb++){
           _p   += we[je]*wnb[je][jnb]*pre[index(ixe+je,ixnb[je]+jnb)] ;
           _T   += we[je]*wnb[je][jnb]*T[index(ixe+je,ixnb[je]+jnb)] ;
           _mub += we[je]*wnb[je][jnb]*mub[index(ixe+je,ixnb[je]+jnb)] ;
	}
	if(_p<0.) _p = 0. ;
        // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}

void EoS1f::eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) {
    cout << " >>>>  no original table avaliable <<<< " << endl;
}


