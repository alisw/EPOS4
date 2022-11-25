//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

class EoS;
class EoSaux;
class EoSChiral : public EoS {
 private:
  EoSaux *eossmall, *eosbig;

 public:
  EoSChiral(char* fileSmall, char* fileLarge);
  ~EoSChiral(void);
  virtual void eosranges(double &_emax, double &_e0, double &_nmax, double &_n0, int &_ne, int &_nn) ;
  virtual void eos(double e, double nb, double nq, double ns, double &_T,
                   double &_mub, double &_muq, double &_mus, double &_p);
  virtual double p(double e, double nb, double nq, double ns);
  virtual void eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p);
};
