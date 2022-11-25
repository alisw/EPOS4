//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#ifndef TRE_H
#define TRE_H

#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;
//==============================================================================
   class Eptree
//==============================================================================
{
public:
  Eptree(){};//default constructor
  //tree,file,#events
  TTree *EposHeadTree;          //epos head tree
  TTree *EposTree;              //epos event tree
  Int_t treEve;                //Number of events in tree
  TFile* fEposTreeFile;

  //header tree
  Int_t iversn;
  Int_t laproj;
  Int_t maproj;
  Int_t latarg;
  Int_t matarg;

  Float_t engy  ;
  Int_t nfreeze;
  Int_t nfull;
  //epos tree
  Int_t      np ; // number of particles                
  Float_t    bim; // centrality variable       
  Float_t *  zus; // different meaning depending on ptl type:
                  // partons: presently unused
                  // hadrons:  decay information :
                  // -999 : hadron is decay product from decay 
                  //        in cascade part (mother unknown)
                  //   -1 : hadron is decay product, mother not stored
                  //   >0 : hadron is decay product, mother index = zus
                  //   -2 : no mother  
  Float_t *  px ; //  -+ 
  Float_t *  py ; //   |  particle
  Float_t *  pz ; //  -+  momentum
  Float_t *  en ; //  energy (version < 3.127) or mass (>= 3.127)
  Float_t *  x  ; //  -+ 
  Float_t *  y  ; //   | particle
  Float_t *  z  ; //   | four position
  Float_t *  t  ; //  -+
  Int_t   *  id ; // particle id (epos code)
  Int_t   *  ist; // status and particle type
                  //  0 and 1 ... hadrons (0 = last generation)
                  //  21 ........ partons
                  //  25 ........ intermediate out-Born partons
                  //  29 ........ string
  Int_t   *  ity; // type of origin
  Int_t   *  ior; // origin (particle index)
  Int_t   *  jor; // origin (particle index)
                  // ior>0, jor=0 : ior is origin
                  // ior>0, jor>0 : origins from  ior to jor
                  // ior=0, jor=0 : no origin
        
        
  Eptree(Int_t ndim){
    //cout<<"        Allocate memory for Eptree arrays "<<endl;
    Xdimension=ndim;
    
    zus = new Float_t [ndim]; 
    px  = new Float_t [ndim];  
    py  = new Float_t [ndim];  
    pz  = new Float_t [ndim];  
    en  = new Float_t [ndim];   
    x   = new Float_t [ndim] ;  
    y   = new Float_t [ndim] ;  
    z   = new Float_t [ndim] ;  
    t   = new Float_t [ndim] ;  
    id  = new Int_t [ndim];  
    ist = new Int_t [ndim]; 
    ity = new Int_t [ndim]; 
    ior = new Int_t [ndim]; 
    jor = new Int_t [ndim]; 
    
  };
  //destructor:
   ~Eptree(){
     //for(Int_t i=0;i<Xdimension;i++){
     //delete id[i];
     // }
     //cout<<"Eptree::~Eptree() "<<px[Xdimension-1]<<" teefile:"<<fEposTreeFile<<endl;
  delete  [] zus;
  delete  [] px ;
  delete  [] py ;
  delete  [] pz ;
  delete  [] en ;
  delete  [] x  ;
  delete  [] y  ;
  delete  [] z  ;
  delete  [] t  ;
  delete  [] id ;
  delete  [] ist;
  delete  [] ity;
  delete  [] ior;
  delete  [] jor;
  fEposTreeFile->TFile::Close();
  //operator delete(zus);
  //cout<<"Eptree::~Eptree() FINISHED"<<endl;
  //cout<<"Delete memory allocated for Eptree arrays"<<endl;
  };
  //function to init tree:
  //virtual void InitTreeFile(char* fnamein, Eptree * eee);
private:
  Int_t Xdimension;
protected:  
#ifdef __ROOT__
  ClassDef(Eptree,1)
#endif
};

void InitTreeFile(char* fnamein, Eptree * eee);

extern "C" void checktreefile_(char* fnamein, int* ierr);

extern "C" void monitorfile_(char *fnamemtr);

//==============================================================================
   class Glob
//==============================================================================
{
public:
        static int Sho;
protected:

#ifdef __ROOT__
  ClassDef(Glob,1)
#endif
};

#endif
