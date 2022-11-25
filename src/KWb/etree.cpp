//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "etree.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>

const int ntotmax = 350000 ;

ETree *etree;

extern "C" void etreeopen_(char* fnamein){
        etree=new ETree(fnamein);
}
extern "C" void etreeclose_(){
        delete etree;
}
extern "C" void etreenevent_(int * nevents){
        *nevents = etree->EposTree->GetEntries();
}
extern "C" void etreeevt_(int * ii, int *_np, float *_bim, float *_sigtot, int *_iextree
           , int *_nev, int *_npt, int *_ngl, int *_kol, int *_nhard, int *_npartproj, int *_nparttarg, int *_nspecp, int *_nspecn //*JJ
           , float *_phi, float *_phir, float *_psi2, float *_psi3, float *_psi4, float *_psi5, float *_ecci2, float *_ecci3, float *_ecci4, float *_ecci5
           , int *_id, int *_ist, int *_ity, int *_ior, int *_jor, float *_px, float *_py, float *_pz, float *_en
           , float *_x, float *_y, float *_z, float *_t) {
        int i=*ii-1;
        int iextree=*_iextree;
        //cout<<" event "<<i<<endl;
    //**JJ std::cout << " ****** AJEDNEK TUATU ******** " << std::endl;
        etree->EposTree->SetBranchAddress("np", _np   );
        etree->EposTree->SetBranchAddress("bim", _bim );
        etree->EposTree->SetBranchAddress("nev", _nev ); //----new---->
        etree->EposTree->SetBranchAddress("npt", _npt );
        etree->EposTree->SetBranchAddress("ngl", _ngl );
        etree->EposTree->SetBranchAddress("kol", _kol );
        if(iextree>=1) //extension 1
        {
         etree->EposTree->SetBranchAddress("sigtot", _sigtot );//*JJ
         etree->EposTree->SetBranchAddress("nhard", _nhard );//*JJ
         etree->EposTree->SetBranchAddress("npartproj", _npartproj );//*JJ
         etree->EposTree->SetBranchAddress("nparttarg", _nparttarg );//*JJ
         etree->EposTree->SetBranchAddress("nspecp", _nspecp );//*JJ
         etree->EposTree->SetBranchAddress("nspecn", _nspecn );//*JJ
        }
        etree->EposTree->SetBranchAddress("phi", _phi );
        etree->EposTree->SetBranchAddress("phir", _phir );
        etree->EposTree->SetBranchAddress("psi2", _psi2 );
        etree->EposTree->SetBranchAddress("psi3", _psi3 );
        etree->EposTree->SetBranchAddress("psi4", _psi4 );
        etree->EposTree->SetBranchAddress("psi5", _psi5 );
        etree->EposTree->SetBranchAddress("ecci2", _ecci2 );
        etree->EposTree->SetBranchAddress("ecci3", _ecci3 );
        etree->EposTree->SetBranchAddress("ecci4", _ecci4 );
        etree->EposTree->SetBranchAddress("ecci5", _ecci5 ); //----------- 
        etree->EposTree->SetBranchAddress("px",  _px  );
        etree->EposTree->SetBranchAddress("py",  _py  );
        etree->EposTree->SetBranchAddress("pz",  _pz  );
        etree->EposTree->SetBranchAddress("e",   _en  );
        etree->EposTree->SetBranchAddress("x",   _x   );
        etree->EposTree->SetBranchAddress("y",   _y   );
        etree->EposTree->SetBranchAddress("z",   _z   );
        etree->EposTree->SetBranchAddress("t",   _t   );
        etree->EposTree->SetBranchAddress("id",  _id  );
        etree->EposTree->SetBranchAddress("ist", _ist );
        etree->EposTree->SetBranchAddress("ity", _ity );
        etree->EposTree->SetBranchAddress("ior", _ior );
        etree->EposTree->SetBranchAddress("jor", _jor );
        etree->EposTree->GetEntry(i);
}
extern "C" void etreehead_(int * _iversn , int * _laproj , int * _maproj , int * _latarg , int * _matarg , float * _engy , int * _nfreeze  ,int * _nfull) {
        etree->EposHeadTree->SetBranchAddress("iversn" , _iversn  );
        etree->EposHeadTree->SetBranchAddress("laproj" , _laproj  );
        etree->EposHeadTree->SetBranchAddress("maproj" , _maproj  );
        etree->EposHeadTree->SetBranchAddress("latarg" , _latarg  );
        etree->EposHeadTree->SetBranchAddress("matarg" , _matarg  );
        etree->EposHeadTree->SetBranchAddress("engy"   , _engy    );
        etree->EposHeadTree->SetBranchAddress("nfreeze", _nfreeze );
        etree->EposHeadTree->SetBranchAddress("nfull"  , _nfull   );
        etree->EposHeadTree->GetEntry(0);
}

