//STARTHEADER
// $Id: fastjetfortran.cc 1570 2009-05-25 10:45:18Z salam $
//
// Copyright (c) 2005-2007, Matteo Cacciari, Gavin Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------
//ENDHEADER

#include <iostream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

namespace fj = fastjet;
using namespace std;

extern "C" {   

// f77 interface to the pp generalised-kt (sequential recombination)
// algorithms, as defined in arXiv.org:0802.1189, which includes
// kt, Cambridge/Aachen and anti-kt as special cases.
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE FASTJETPPSEQREC(P,NPART,R,PALG,F77JETS,NJETS)
//   DOUBLE PRECISION P(4,*), R, PALG, F, F77JETS(4,*)
//   INTEGER          NPART, NJETS
// 
// where on input
//
//   P        the input particle 4-momenta
//   NPART    the number of input momenta
//   R        the radius parameter
//   PALG     the power for the generalised kt alg 
//            (1.0=kt, 0.0=C/A,  -1.0 = anti-kt)
//
// and on output 
//
//   F77JETS  the output jet momenta (whose second dim should be >= NPART)
//            sorted in order of decreasing p_t.
//   NJETS    the number of output jets 
//
// For the values of PALG that correspond to "standard" cases (1.0=kt,
// 0.0=C/A, -1.0 = anti-kt) this routine actually calls the direct
// implementation of those algorithms, whereas for other values of
// PALG it calls the generalised kt implementation.
//
void fastjetppgenkt_(const double * p, const int & npart,                   
                     const double & R, const double & palg,
                     double * f77jets, int & njets,
                     int * array1, int * array2,  const int & ipri) {

if(ipri)cout <<"initialize fastjetppgenkt_ "<<endl;
    
    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      fj::PseudoJet psjet(mom);
      psjet.set_user_index(i+1) ; 
      input_particles.push_back(psjet);    
    }
    
    // prepare jet def and run fastjet
    fj::JetDefinition jet_def;
    if (palg == 1.0) {
      jet_def = fj::JetDefinition(fj::kt_algorithm, R);
    }  else if (palg == 0.0) {
      jet_def = fj::JetDefinition(fj::cambridge_algorithm, R);
    }  else if (palg == -1.0) {
      jet_def = fj::JetDefinition(fj::antikt_algorithm, R);
    } else {
      jet_def = fj::JetDefinition(fj::genkt_algorithm, R, palg);
    }

if(ipri)cout <<"perform clustering     "<<endl;
    
    fj::ClusterSequence cs(input_particles, jet_def);

if(ipri)cout <<"extract jets  (pt-ordered)    "<<endl;

    vector<fj::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    njets = jets.size();

if(ipri)cout <<"transfer jets to f77jets"<<endl;

    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

if(ipri)cout <<"transfer constits to array1 array2"<<endl;

    int count=0;
    //if(ipri)cout <<"nconstits for the "<<njets<<" jets:"<<endl;
    for (int i=0; i<njets; i++) {
        vector<fj::PseudoJet>  constits  ;
        constits = cs.constituents(jets[i])  ;
        int nconstits  =   constits.size() ;
        //if(ipri)cout<< nconstits<<" ";
        for (int j=0; j<nconstits; j++,count++) {
          int indx = constits[j].user_index() ;
          *array2 = indx;
          array2++;
        }
        *array1=count;
        array1++;
    }
    //if(ipri)cout<<endl;
    // clean up
    
}

}

