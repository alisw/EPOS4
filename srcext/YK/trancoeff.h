//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

class EoS ;

class TransportCoeff
{
 double etaS, zetaS, taupi, tauPi ;
 EoS *eos ;
 public:
 TransportCoeff(double _etaS, double _zetaS, EoS *_eos) ;
 ~TransportCoeff() {} ;
 void getEta(double e, double T, double &_etaS, double &_zetaS) ;
 void getTau(double T, double &_taupi, double &_tauPi) ;
 inline bool isViscous() { if(etaS>0. || zetaS>0.) return true ; else return false ; }
} ;
