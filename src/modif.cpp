//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//


#include <iostream>
extern "C" {
  void getnhpom_(int*);
  void getrng1_(int*);
  void pfe3modif_();
  void pfe3paramset5_(float* v1,float* v2,float* v3,float* v4,float* v5);
}


//------------------------------------------------------------
void pfe3modif_(){
  //------------------------------------------------------------
  
  // Finetuning for a particular system.
  // Uncomment and modify the definitions of key variables.
  // All parameters may (in principle) be made N or Z dependent
  
  // Get N and Z values
  // ------------------

  // int N, Z;
  // float yrmax ,yco, ycoj, fec, taufo;
  // getnhpom_(&N);    // Number of hard Pomerons. Used for pp.   
  // getrng1_(&Z);     // Centrality (0 = 'peripheral' and 1 = 'central'). Used for AA
                        // In rare cases, Z may be slighly bigger than 1
  
  // pp 7TeV
  // -------
  
  //      yrmax=1+0.04*min(N,12)  !radial boost
  //      ycoi= 1.0                !longitudinal boost
  //      ycoj=0.5                 !longitudinal boost rapidity dependence
  //      fecc=min(0.015*N,0.15)  !eccentricity
  //      taufo=1.0               !duration of expansion (in fm/c)
  
  // PbPb 5TeV
  // ---------
  
  //      yrmax=1+0.16*Z                     !radial boost
  //      ycoi= 1.0                          !longitudinal boost
  //      ycoj=0.5                           !longitudinal boost rapidity dependence
  //      fecc=max(0.05,min(0.23,.4-.4*Z))   !eccentricity
  //      taufo=1+5*Z                        !duration of expansion (in fm/c)
  
  // AuAu 39GeV
  // ----------
  
  //      yrmax= 0.80+0.05*Z                  !radial boost
  //      ycoi=1.0                            !longitudinal boost
  //      ycoj=0.5                            !longitudinal boost rapidity dependence
  //      fecc=max(0.065,min(0.23,.29-.26*Z)) !eccentricity
  //      taufo=1+5*Z                         !duration of expansion (in fm/c)
  
  // Apply the values
  //-----------------
  // pfe3paramset5_(&yrmax ,&ycoi, &ycoj, &fecc, &taufo);

}     
