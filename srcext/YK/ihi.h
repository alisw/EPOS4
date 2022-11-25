//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

extern "C" void open2dhisto_(char * filename, int * kmax, int * nxbins, int * nybins, double * etamax) ;
extern "C" void fill2dhisto1_(int * kk, double * delta_eta, double * delta_phi) ;
extern "C" void fill2dhisto2_(int * kk, double * delta_eta, double * delta_phi) ;
extern "C" void close2dhisto_() ;
