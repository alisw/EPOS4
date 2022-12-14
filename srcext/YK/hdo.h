//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//


class Cell ;
class Fluid ;
class EoS ;
class TransportCoeff ;

class Hydro
{
private :
	Fluid *f ;
	EoS *eos ;
	TransportCoeff *trcoeff ;
	double dt, tau ;
	double tau_z ;
public :
	Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0, double _dt) ;
	~Hydro();
	void setDtau(double deltaTau) ;
 void setFluid(Fluid *_f) { f = _f; }
 Fluid* getFluid() { return f; }

	void hlle_flux(Cell *left, Cell *right, int direction, int mode) ;
	void visc_flux(Cell *left, Cell *right, int direction) ;
	void visc_source_step(int ix, int iy, int iz) ;
	void source(double tau, double x, double y, double z, double Q[7], double S[7]) ;
	void source_step(int ix, int iy, int iz, int mode) ;
	void NSquant(int ix, int iy, int iz, double pi[][4], double& Pi, double dmu[4][4], double &du) ;
	void setNSvalues() ;
	void ISformal() ;
	void setQfull() ;
	void performStep(void) ;
	inline double getTau() { return tau ; }
};
