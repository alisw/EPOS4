//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#pragma once
#include <math.h>
#include "inc.h"
class TGraph ;

class EoS
{
public :
	virtual ~EoS() { return ; }
	virtual void eosranges(double &emax, double &e0, double &nmax, double &n0, int &ne, int &nn) = 0 ;
	virtual void eos(double e, double nb, double nq, double ns,
		double &T, double &mub, double &muq, double &mus, double &p) = 0 ;
	virtual double p(double e, double nb, double ns, double nq) = 0;
    virtual void eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) = 0;
    double s(double e, double nb, double nq, double ns) ;
	inline double cs2(void) { return 1./3. ; }
	inline double cs(void) { return sqrt(1./3.) ; }
	virtual inline double cs2(double e) { return 1./3. ; };

};

class EoSs : public EoS
{
private:
 TGraph *gp, *gT, *gmu ;
public:
	EoSs(string fname, int ncols);
	~EoSs();

	virtual inline void eosranges(double &emax, double &e0, double &nmax, double &n0, int &ne, int &nn) {}
	virtual inline void eos(double e, double nb, double nq, double ns,
		double &T, double &mub, double &muq, double &mus, double &_p){
			_p = p(e) ;
			T = t(e) ;
			mub = muq = mus = 0. ;
	}
	virtual inline double p(double e, double nb, double ns, double nq)
	{ return p(e) ; }

	double p(double e) ;
	double dpe(double e) ;
	double t(double e) ;
	double mu(double e) ;
    virtual void eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) =0;
	virtual double cs2(double e) { return dpe(e) ; }
	//virtual double cs(double e) { return sqrt(dpe(e)) ; }

};
