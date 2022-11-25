//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "rmn.h"
#include <iomanip>

extern void handler(int sig) ;

bool debugRiemann ;

/*
void transformPVold(EoS *eos, double Q[7], double &e, double &p, double &nb, double &nq, double &ns, double &vx, double &vy, double &vz)
{
	const double dpe = 1./3. ;
	const double corrf = 0.9999 ;
	double M = sqrt(Q[X_]*Q[X_]+Q[Y_]*Q[Y_]+Q[Z_]*Q[Z_]) ;
	double v = 0.5, v0 = 0.5 ;
	if(M==0.){
		e=Q[T_]; vx=0.; vy=0.; vz=0.; nb=Q[NB_]; nq=Q[NQ_]; ns=Q[NS_];
		p=eos->p(e,nb,nq,ns) ;
		return ;
	}
	if(Q[T_]<0.) { e=0.; p = 0.; vx = vy = vz = 0.; nb=nq=ns=0.; return ; }
	if(M>Q[T_]){
		Q[X_] *= corrf*Q[T_]/M ;
		Q[Y_] *= corrf*Q[T_]/M ;
		Q[Z_] *= corrf*Q[T_]/M ;
		M = Q[T_]*corrf ;
	}
	int i = 0 ;  //debug
	do{
		v0 = v ;
		e = Q[T_]-M*v0 ;
		if(e<0.) e = 0. ;
		nb = Q[NB_]*sqrt(1-v*v) ;
		nq = Q[NQ_]*sqrt(1-v*v) ;
		ns = Q[NS_]*sqrt(1-v*v) ;
		p = eos->p(e, nb, nq, ns) ;
		v = v0 - ( M - (Q[T_] + p)*v0 ) / ( -(Q[T_]+p) + M*v0*dpe ) ;
		if(v>1.) v=0.99999 ;
		if(v<0.) v = 0. ;
		i++ ;
//		if(i>1000 && i<1005) cout << "inf.loop in EoS: e= " << setw(14) << e << " nb= " << setw(12) << nb
//			<< " nq= " << setw(12) << nq << " ns= " << setw(12) << ns << endl ;
		if(i>1010) break ;
	}
	while(fabs(v-v0)>0.000001) ;
	// alternative procedure
	if(i>1000){
		double vmin=0., vmax=1. ;
		do{
			v = 0.5*(vmin+vmax) ;
			e = Q[T_] - M*v ;
			nb = Q[NB_]*sqrt(1-v*v) ;
			nq = Q[NQ_]*sqrt(1-v*v) ;
			ns = Q[NS_]*sqrt(1-v*v) ;
			p = eos->p(e, nb, nq, ns) ;
			if ( M/(Q[T_]+p)-v > 0.) vmin = v ;
			else vmax = v ;
			}while((vmax-vmin)>0.00001) ;
//		cout << "solved in EoS: e= " << setw(12) << e << " nb= " << setw(12) << nb
//			<< " nq= " << setw(12) << nq << " ns= " << setw(12) << ns << endl ;
	}
	vx = v*Q[X_]/M ;
	vy = v*Q[Y_]/M ;
	vz = v*Q[Z_]/M ;
	e = Q[T_] - M*v ;
	p = eos->p(e,nb,nq,ns) ;
	nb = Q[NB_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	nq = Q[NQ_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	ns = Q[NS_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	if(e<0. || sqrt(vx*vx+vy*vy+vz*vz)>1.){
		cout << Q[T_] << "  " << Q[X_] << "  " << Q[Y_] << "  " << Q[Z_] << "  " << Q[NB_] << endl ;
		cout<<"transformRF::Error" ;
	}
	if(!(nb<0. || nb>=0.)){
		cout<<"transformRF::Error nb=#ind" ;
		return ;
	}
}
*/


void transformPV(EoS *eos, double Q[7], double &e, double &p, double &nb, double &nq, double &ns, double &vx, double &vy, double &vz)
{
	const int MAXIT = 100 ;
	const double dpe = 1./3. ;
	const double corrf = 0.9999 ;
	double v, vl = 0., vh = 1., dvold, dv, f, df ;
	if(debugRiemann){
	cout << "Riemann debug---------------\n" ;
	cout << setw(14) << Q[0] << setw(14) << Q[1] << setw(14) << Q[2] << setw(14) << Q[3] << endl ;
	cout << setw(14) << Q[4] << setw(14) << Q[5] << setw(14) << Q[6] << endl ;
	}
	double M = sqrt(Q[X_]*Q[X_]+Q[Y_]*Q[Y_]+Q[Z_]*Q[Z_]) ;
	if(Q[T_]<=0.) { e=0.; p = 0.; vx = vy = vz = 0.; nb=nq=ns=0.; return ; }
	if(M==0.){
		e=Q[T_]; vx=0.; vy=0.; vz=0.; nb=Q[NB_]; nq=Q[NQ_]; ns=Q[NS_];
		p=eos->p(e,nb,nq,ns) ; if(p==-999.)cout<<"WARNING 1 in transformPV: p = "<<p<<" (out of range)  e nb nq ns = "<<e<<" "<<nb<<" "<<nq<<" "<<ns<<endl;
		return ;
	}
	if(M>Q[T_]){
		Q[X_] *= corrf*Q[T_]/M ;
		Q[Y_] *= corrf*Q[T_]/M ;
		Q[Z_] *= corrf*Q[T_]/M ;
		M = Q[T_]*corrf ;
	}

	v = 0.5*(vl+vh) ; 
	e = Q[T_]-M*v ;
    if(e<0.) e = 0. ;
	nb = Q[NB_]*sqrt(1-v*v) ;
	nq = Q[NQ_]*sqrt(1-v*v) ;
	ns = Q[NS_]*sqrt(1-v*v) ;
	p = eos->p(e, nb, nq, ns) ; if(p==-999.)cout<<"WARNING 2 in transformPV: p = "<<p<<" (out of range)  e nb nq ns = "<<e<<" "<<nb<<" "<<nq<<" "<<ns<<endl;
	f  = (Q[T_] + p) * v - M   ;
	df = (Q[T_] + p) - M*v*dpe ;
	dvold = vh -vl ;
	dv = dvold ;
	for(int i=0; i<MAXIT; i++){
		if( (f+df*(vh-v))*(f+df*(vl-v)) > 0. || fabs(2.*f) > fabs(dvold*df) ){ // bisection
			dvold = dv ;
			dv = 0.5*(vh-vl) ;
			v = vl + dv ;
//			cout << "BISECTION v = " << setw(12) << v << endl ;
		}else{ // Newton
			dvold = dv ;
			dv = f / df ;
			v -= dv ;
//			cout << "NEWTON v = " << setw(12) << v << endl ;
		}
		if(fabs(dv)<0.00001) break ;

		e = Q[T_]-M*v ;
		if(e<0.) e = 0. ;
     //   if (v*v >= 1){
     //       cout << " ********************** " << endl;
     //       cout << " HERE WE HAVE A PROBLEM v= " << v << "for i: " <<i << "  f:  " << f << "   df:   " << df << endl;
     //       cout << " p= " << p << " Q[T_]= " << Q[T_] << " T_= " << T_ << " M= " << M << " e= " << e << " nb= " << nb << endl;
     //       cout << " ********************** " << endl;
     //   }
		nb = Q[NB_]*sqrt(1-v*v) ;
		nq = Q[NQ_]*sqrt(1-v*v) ;
		ns = Q[NS_]*sqrt(1-v*v) ;
		p = eos->p(e, nb, nq, ns) ; if(p==-999.)cout<<"WARNING 3 in transformPV: p = "<<p<<" (out of range)  e nb nq ns = "<<e<<" "<<nb<<" "<<nq<<" "<<ns<<endl;
		f = (Q[T_] + p)*v - M ;
		df = (Q[T_] + p) - M*v*dpe ;

		if(f>0.) vh = v ;
		else vl = v ;
		if(debugRiemann)
		 cout << "step " << i << "  " << e << "  " << nb << "  " << nq << "  " << ns << "  " << p
		  << "  " << vx << "  " << vy << "  " << vz << endl ;
//		if(i==40) { cout << "error : does not converge\n" ; exit(1) ; } ;
	} // for loop
//----------after
//	v = 0.5*(vh+vl) ;
	vx = v*Q[X_]/M ;
	vy = v*Q[Y_]/M ;
	vz = v*Q[Z_]/M ;
	e = Q[T_] - M*v ;
	p = eos->p(e,nb,nq,ns) ;   if(p==-999.)cout<<"WARNING 4 in transformPV: p = "<<p<<" (out of range)  e nb nq ns = "<<e<<" "<<nb<<" "<<nq<<" "<<ns<<endl;
	nb = Q[NB_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	nq = Q[NQ_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	ns = Q[NS_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	if(e<0. || sqrt(vx*vx+vy*vy+vz*vz)>1.){
		cout << Q[T_] << "  " << Q[X_] << "  " << Q[Y_] << "  " << Q[Z_] << "  " << Q[NB_] << endl ;
		cout<<"transformRF::Error\n" ;
	}
	if(!(nb<0. || nb>=0.)){
		cout<<"transformRF::Error nb=#ind\n" ;
//		return ;
	}
}


void transformPVBulk(EoS *eos, double Pi, double Q[7], double &e, double &p, double &nb, double &nq, double &ns, double &vx, double &vy, double &vz)
{
	const int MAXIT = 100 ;
	const double dpe = 1./3. ;
	const double corrf = 0.9999 ;
	double v, vl = 0., vh = 1., dvold, dv, f, df ;
	if(debugRiemann){
	cout << "Riemann debug---------------\n" ;
	cout << setw(14) << Q[0] << setw(14) << Q[1] << setw(14) << Q[2] << setw(14) << Q[3] << endl ;
	cout << setw(14) << Q[4] << setw(14) << Q[5] << setw(14) << Q[6] << endl ;
	}
	double M = sqrt(Q[X_]*Q[X_]+Q[Y_]*Q[Y_]+Q[Z_]*Q[Z_]) ;
	if(Q[T_]<=0.) { e=0.; p = 0.; vx = vy = vz = 0.; nb=nq=ns=0.; return ; }
	if(M==0.){
		e=Q[T_]; vx=0.; vy=0.; vz=0.; nb=Q[NB_]; nq=Q[NQ_]; ns=Q[NS_];
		p=eos->p(e,nb,nq,ns)+Pi ;
		return ;
	}
	if(M>Q[T_]){
		Q[X_] *= corrf*Q[T_]/M ;
		Q[Y_] *= corrf*Q[T_]/M ;
		Q[Z_] *= corrf*Q[T_]/M ;
		M = Q[T_]*corrf ;
	}

	v = 0.5*(vl+vh) ;
	e = Q[T_]-M*v ;
	if(e<0.) e = 0. ;
	nb = Q[NB_]*sqrt(1-v*v) ;
	nq = Q[NQ_]*sqrt(1-v*v) ;
	ns = Q[NS_]*sqrt(1-v*v) ;
	p = eos->p(e, nb, nq, ns)+Pi ;
	f = (Q[T_] + p)*v - M ;
	df = (Q[T_] + p) - M*v*dpe ;
	dvold = vh -vl ;
	dv = dvold ;
	for(int i=0; i<MAXIT; i++){
		if( (f+df*(vh-v))*(f+df*(vl-v)) > 0. || fabs(2.*f) > fabs(dvold*df) ){ // bisection
			dvold = dv ;
			dv = 0.5*(vh-vl) ;
			v = vl + dv ;
//			cout << "BISECTION v = " << setw(12) << v << endl ;
		}else{ // Newton
			dvold = dv ;
			dv = f / df ;
			v -= dv ;
//			cout << "NEWTON v = " << setw(12) << v << endl ;
		}
		if(fabs(dv)<0.00001) break ;

		e = Q[T_]-M*v ;
		if(e<0.) e = 0. ;
		nb = Q[NB_]*sqrt(1-v*v) ;
		nq = Q[NQ_]*sqrt(1-v*v) ;
		ns = Q[NS_]*sqrt(1-v*v) ;
		p = eos->p(e, nb, nq, ns)+Pi ;
		f = (Q[T_] + p)*v - M ;
		df = (Q[T_] + p) - M*v*dpe ;

		if(f>0.) vh = v ;
		else vl = v ;
		if(debugRiemann)
		 cout << "step " << i << "  " << e << "  " << nb << "  " << nq << "  " << ns << "  " << p
		  << "  " << vx << "  " << vy << "  " << vz << endl ;
//		if(i==40) { cout << "error : does not converge\n" ; exit(1) ; } ;
	} // for loop
//----------after
//	v = 0.5*(vh+vl) ;
	vx = v*Q[X_]/M ;
	vy = v*Q[Y_]/M ;
	vz = v*Q[Z_]/M ;
	e = Q[T_] - M*v ;
	p = eos->p(e,nb,nq,ns) ;
	nb = Q[NB_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	nq = Q[NQ_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	ns = Q[NS_]*sqrt(1-vx*vx-vy*vy-vz*vz) ;
	if(e<0. || sqrt(vx*vx+vy*vy+vz*vz)>1.){
		cout << Q[T_] << "  " << Q[X_] << "  " << Q[Y_] << "  " << Q[Z_] << "  " << Q[NB_] << endl ;
		cout<<"transformRF::Error\n" ;
	}
	if(!(nb<0. || nb>=0.)){
		cout<<"transformRF::Error nb=#ind\n" ;
		cout<<"Q [7]: "<<Q[0]<<" "<<Q[1]<<" "<<Q[2]<<" "<<Q[3]<<" "<<Q[4]<<" "<<Q[5]<<" "<<Q[6]<<endl ;
		cout<<"e vx vy vz nb nq ns: "<<e<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<nb<<" "<<nq<<" "<<ns<<endl ;
		handler(333) ;
//		return ;
	}
}


void transformCV(double e, double p, double nb, double nq, double ns, double vx, double vy, double vz, double Q [])
{
	const double gamma2 = 1./(1-vx*vx-vy*vy-vz*vz) ;
	Q[T_] = (e+p*(vx*vx+vy*vy+vz*vz))*gamma2 ;
	Q[X_] = (e+p)*vx*gamma2 ;
	Q[Y_] = (e+p)*vy*gamma2 ;
	Q[Z_] = (e+p)*vz*gamma2 ;
	Q[NB_] = nb*sqrt(gamma2) ;
	Q[NQ_] = nq*sqrt(gamma2) ;
	Q[NS_] = ns*sqrt(gamma2) ;
}
