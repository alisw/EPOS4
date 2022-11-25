//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "eos.h"

class EoS3f : public EoS
{
private:
 double emax, nmax, e0, n0 ;
 int ne, nn ;
 double B , volex0 , delta0 , aaa , bbb;
 double *egrid, **ngrid ;
 double *T, *pre, *mub, *muq, *mus ;
 inline int index(int ie, int ib, int iq, int is){ return ie + ne*ib + ne*nn*iq + ne*nn*nn*is ; }
public:
 EoS3f(char *filename, double B, double volex0, double delta0, double aaa, double bbb);
 ~EoS3f(void);

 virtual void eosranges(double &_emax, double &_e0, double &_nmax, double &_n0, int &_ne, int &_nn) ;
 void getue(double e, int &ixe, double &ue, int &iout) ; 
 void getun(int ie, double n, int &ixn_low, int &ixn_up,
                    double &un_low, double &un_up, int &iout) ;
 virtual void eos(double e, double nb, double nq, double ns,
         double &_T, double &_mub, double &_muq, double &_mus, double &_p) ;
 virtual inline double p(double e, double nb, double nq, double ns)
 { double T, mub, muq, mus, pp ; eos(e, nb, nq, ns, T, mub, muq, mus, pp) ; return pp ; }
  virtual void eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p);
 friend void mix(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p);
};
