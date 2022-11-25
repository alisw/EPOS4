//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "trancoeff.h"
#include "eos.h"

TransportCoeff::TransportCoeff(double _etaS, double _zetaS, EoS *_eos)
{
 etaS = _etaS ;
 zetaS = _zetaS ;
 eos = _eos ;
}

void TransportCoeff::getEta(double e, double T, double &_etaS, double &_zetaS)
{
 _etaS = etaS ;
 _zetaS=zetaS*(1./3.-eos->cs2(e))/(exp((0.16-T)/0.001)+1.) ;
}


void TransportCoeff::getTau(double T, double &_taupi, double &_tauPi)
{
	if(T>0.) _taupi=max(3./5.068*etaS/T,0.3) ; else _taupi=0.1 ;
	if(T>0.) _tauPi=max(3./5.068*(1./4./C_PI)/T,0.05) ; else _tauPi=0.1 ;
}
