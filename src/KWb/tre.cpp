//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <TTree.h>
#include <TFile.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "tre.h"
  
using namespace std;

//==============================================================================
//void Eptree::InitTreeFile(char* fnamein, Eptree * eee)
void InitTreeFile(char* fnamein, Eptree * eee)
//==============================================================================
{
  TTree *hd;
  TTree *ep;
  int fSho=0; //Glob::Sho;
  if(fSho)cout<<"---------------- Init Epos Trees ----------------------"<<endl;
  //if(fSho)cout<<"Open tree file " << fnamein <<"  ";
  eee->fEposTreeFile = TFile::Open(fnamein);
  //cout<<"file "<<eee->fEposTreeFile<<endl;
  if(fSho)cout <<endl;
  if(eee->fEposTreeFile==NULL) {
    cout<<"Error: can't open file "<< fnamein <<endl;
    exit(-1);
  } 
  hd=(TTree*)eee->fEposTreeFile->Get("teposhead");
  eee->EposHeadTree=(TTree*)hd;
  hd->SetBranchAddress("iversn", &eee->iversn);
  hd->SetBranchAddress("laproj", &eee->laproj);
  hd->SetBranchAddress("maproj", &eee->maproj);
  hd->SetBranchAddress("latarg", &eee->latarg);
  hd->SetBranchAddress("matarg", &eee->matarg);
  hd->SetBranchAddress("engy",   &eee->engy  );         
  hd->SetBranchAddress("nfreeze",&eee->nfreeze);    
  hd->SetBranchAddress("nfull",  &eee->nfull ); 
  if(fSho)cout<<"hd tree entries : "<<hd->GetEntries()<<endl;
  hd->GetEntry(0);
  if(fSho)cout<<"version                : "<< eee->iversn <<endl;
  if(fSho)cout<<"laproj                 : "<< eee->laproj <<endl;
  if(fSho)cout<<"maproj                 : "<< eee->maproj <<endl;
  if(fSho)cout<<"latarg                 : "<< eee->latarg <<endl;
  if(fSho)cout<<"matarg                 : "<< eee->matarg <<endl;
  if(fSho)cout<<"engy                   : "<< eee->engy   <<endl;
  if(fSho)cout<<"nfreeze                : "<< eee->nfreeze<<endl;
  if(fSho)cout<<"nfull                  : "<< eee->nfull  <<endl;
  ep=(TTree*)eee->fEposTreeFile->Get("teposevent");
  eee->EposTree=(TTree*)ep;
  eee->treEve=ep->GetEntries();
  if(fSho)cout<<"ep tree entries : "<<eee->treEve<<endl;
  //cout << eee->px << endl;
  ep->SetBranchAddress("np",  &eee->np);
  ep->SetBranchAddress("bim", &eee->bim);
  ep->SetBranchAddress("zus", eee->zus);
  ep->SetBranchAddress("px",  eee->px);
  ep->SetBranchAddress("py",  eee->py);
  ep->SetBranchAddress("pz",  eee->pz);
  ep->SetBranchAddress("e",   eee->en);
  ep->SetBranchAddress("x",   eee->x);
  ep->SetBranchAddress("y",   eee->y);
  ep->SetBranchAddress("z",   eee->z);
  ep->SetBranchAddress("t",   eee->t);
  ep->SetBranchAddress("id",  eee->id);
  ep->SetBranchAddress("ist", eee->ist);
  ep->SetBranchAddress("ity", eee->ity);
  ep->SetBranchAddress("ior", eee->ior);
  ep->SetBranchAddress("jor", eee->jor);
  //cout << eee->px << endl;
  //ep->GetEntry(0);
  //cout << eee->px << "	" << eee->py << endl;
  if(fSho)cout<<"---------------- End of Init Epos Trees -----------------"<<endl;

}

//======================================================================================
void checktreefile_(char* fnamein, Int_t* ierr) 
//======================================================================================
{
  //KM: fTreeFile to tTreeFile, fEposHeadTree to tEposHeadTree, fEposTree to tEposTree
  //Please, do not use global variable in local tests
  TFile *tTreeFile;
  TTree* tEposHeadTree; 
  TTree* tEposTree;   
  tTreeFile = TFile::Open(fnamein);
  if(tTreeFile==NULL) {
    cout<<"Skip "<< fnamein << " (No file)" <<  endl;
    *ierr=1;
    return;  
  }
  if(tTreeFile->IsZombie()){
    cout<<"Skip "<< fnamein << " (Zombi)" <<  endl;
    *ierr=1;
    tTreeFile->Close(); return;  
  }
  TTree* tree ;
  tTreeFile->GetObject("teposhead;1",tree) ;
  if(!tree){
    cout<<"Skip "<< fnamein << " (No tree teposhead)" <<  endl;
    *ierr=1;
    tTreeFile->Close(); return;  
  }
  tree->Delete();
  TTree* tree2 ;
  tTreeFile->GetObject("teposevent;1",tree2) ;
  if(!tree2){
    cout<<"Skip "<< fnamein << " (No tree teposevent)" <<  endl;
    *ierr=1;
    tTreeFile->Close(); return;  
  }
  tree2->Delete();
  
  tEposHeadTree=(TTree*)tTreeFile->Get("teposhead");
  tEposTree=(TTree*)tTreeFile->Get("teposevent");
  if(tEposHeadTree == NULL) {
    cout<<"Sip "<<fnamein << "(teposhead does not exist)" <<  endl;
    *ierr=1;
    tTreeFile->Close(); return;  
  }
  if(tEposTree == NULL) {
    cout<<"Sip "<<fnamein << "(teposevent does not exist)" <<endl;
    *ierr=1;
    tTreeFile->Close(); return;  
  }  
  *ierr=0; 
  tTreeFile->Close();  
}   

//======================================================================================
void monitorfile_(char *ffnamemtr)
//======================================================================================
{
  ofstream mfile ;
  if(strlen(ffnamemtr)>0) {
    mfile.open(ffnamemtr,ios_base::app) ;
    cout<<"ffnamemtr opened"<<endl;
    cout.rdbuf(mfile.rdbuf()) ;
    cerr.rdbuf(mfile.rdbuf()) ;
    cout<<"ffnamemtr is new stdout, stderr"<<endl;
   }
}
 


//Eptree::~Eptree(){
  //for(Int_t i=0;i<Xdimension;i++){
  // delete zus[i];
  //}
//  cout<<"Eptree::~Eptree() do nothing "<<t[10]<<endl;
  //delete  *zus[0];
  //operator delete(zus);
//}
