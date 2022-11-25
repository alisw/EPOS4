//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

extern "C" void eosihlle_(double *T, double *mu_b, double *mu_q, double *mu_s, double *e, double *n_b, double *n_q, double *n_s, double *p) ;
extern "C" void eosohlle_(double *e, double *n_b, double *n_q, double *n_s, double *T, double *mu_b, double *mu_q, double *mu_s, double *p) ;
extern "C" void loadptldky_(char *fileparticles, char *filedecay) ;
extern "C" void checkdecays_(void) ;
extern "C" void setparameters_(double *_gammaS, double *_B, double *_mS, double *_volex0, double *_delta0, double *_aaa, double *_bbb, double *_muc, double *_emax, double *_e0,
	double *_anmax, double *_bnmax, double *_an0, int *_nxe, int *_nxn, double *_muBmax, double *_muQmax, double *_muSmax, double *_muBmx, double *_muQmx, double *_muSmx, double *_temmax, int *_maxIter, double *_Taccuracy) ;
extern "C" void seteostype_(int *type) ;
extern "C" void makeinvtable_(char *fnamemtr, char *filename, int *istart, int *iend, int *jstart, int *jend) ;
extern "C" void makeinvtable1f_(char *fnamemtr, char *filename, int *istart, int *iend, int *jstart, int *jend) ;
extern "C" void maketctable_(char *filename, int *istart, int *iend) ;

extern "C" void initeos3tc_(char *filename) ;

extern "C" int igetnparticles_(void) ;
extern "C" void igetparticle_(char* name, double* mass, int* charge, int* id_PDG, 
	int* degeneracy, int* stats ,int* baryon_number, int* strangeness, int* flast_particle, int* flag) ;
extern "C" void destroyptldky_(void) ;

void yields(double T, double mu_b, double mu_q, double mu_s, double volume) ;
void yields_exv2(double T, double mu_b, double mu_q, double mu_s, double volume) ;
void HG_exv2(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) ;

double getT_c(double mu_b, double mu_q, double mu_s) ;

#define EOS3F 22
#define EOSFO 2

class TC
{
private:
        double aa ;
        double bb ; 
	int nn ;
	double *TcTab ;
	inline int index(int ib, int iq, int is){ return ib + nn*iq + nn*nn*is ; }
public:
	TC(char *filename);
	~TC(void);
        void getu(double mu, double muMax, int &i, double &u) ;
        double pi(double mub, double muq, double mus) ;
};
