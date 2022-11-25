//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <TRandom3.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <cstdlib> 
int *id, *ist, *ity, *ior, *jor, *np, *iversn, *laproj, *maproj, *latarg, *matarg, *nfull, *nfreeze ;
float *bim, *engy, *sigtot ; //*JJ
float *zus, *px, *py, *pz, *e, *x, *y, *z, *t;
int *nev, *npt, *ngl, *kol; 
int *nhard, *npartproj, *nparttarg, *nspecp, *nspecn; //*JJ
float *phi, *phir, *psi2, *psi3, *psi4, *psi5, *ecci2, *ecci3, *ecci4, *ecci5;
const int Mxpom = 40 ;
TTree* thead[Mxpom];         
TTree *tree[Mxpom] ;
TFile *file[Mxpom] ;
const int ntotmax = 350000 ;

extern "C" void opentree_(int *_npom, int *_iextree, char* filename){  
//*JJ    std:: cout << "!!!!!!!!!!!!!!!!!!! OTWIERAM DRZEWO !!!!!!!!!!!!!!!!!" << std::endl;
        if(*_npom>Mxpom-1)
        {
           std::cout<<"Error: npom = "<<*_npom<<" is too big"<<std::endl;
           exit(-1);
        } 

        int npom=*_npom-1;
        int iextree=*_iextree;

        char name [200] ; 
        strcpy(name,filename);
        //printf(" opentree *************  %s *************\n",name)  ;
        //std::cout<<" *****fil****  *_npom = "<<npom+1<<std::endl;
        file[npom] = new TFile(filename, "RECREATE"); //output file
        //std::cout << " *****fil**** open tree; file: " << filename << std::endl ;
    
        iversn = new int ;    //EPOS version number
        laproj = new int ;    //atomic number projectile
        maproj = new int ;    //mass number projectile
        latarg = new int ;    //atomic number target
        matarg = new int ;    //mass number target
        engy   = new float ;  //energy in the cms in GeV
        nfull  = new int ;    //number of full events
        nfreeze= new int ;    //number of freeze outs per full event (to increase stat)

        sprintf(name,"teposhead%i",npom);
        //std::cout<<" *****fil**** thead name  "<<name<<std::endl;
        thead[npom] = new TTree(name,"header"); 
        thead[npom]->Branch("iversn",iversn,"iversn/I");
        thead[npom]->Branch("laproj",laproj,"laproj/I");
        thead[npom]->Branch("maproj",maproj,"maproj/I");
        thead[npom]->Branch("latarg",latarg,"latarg/I");
        thead[npom]->Branch("matarg",matarg,"matarg/I");
        thead[npom]->Branch("engy",engy,"engy/F");
        thead[npom]->Branch("nfull",nfull,"nfull/I");
        thead[npom]->Branch("nfreeze",nfreeze,"nfreeze/I");
      
        zus = new float [ntotmax]; //private use
        px = new float [ntotmax];  //  p_x of particle
        py = new float [ntotmax];  //  p_y of particle
        pz = new float [ntotmax];  //  p_z of particle
        e = new float [ntotmax];   //  energy of particle
        x = new float [ntotmax];   //  x component of formation point
        y = new float [ntotmax];   //  y component of formation point
        z = new float [ntotmax];   //  z component of formation point
        t = new float [ntotmax];   //  formation time
        id = new int [ntotmax];    //  particle id  (see array "idt" in function "idtrafo" in file "ids.f": first column "epos id" second one PDG id
        ist = new int [ntotmax];   //  particle status (hadron last generation (0) or not (1); other numbers refer to partons, Pomerons, etc)  
        ity = new int [ntotmax];   //  type of particle origin (20-29 from soft strings, 30-39 from hard strings, 40-59 from remnants, 60 from fluid)
        ior = new int [ntotmax];   //  index of father  (resonance decay products have only a father)
        jor = new int [ntotmax];   //  index of mother  (mothers are needed for exemple for strings: the partons between ior and jor constitute the string)
        np = new int ;          //  number of particles
        bim = new float ;       //  impact parameter (usually; other choices are possible)
        nev = new int ;  //----new---->
        npt = new int ; 
        ngl = new int ; 
        kol = new int ;
        if(iextree>=1) //extension 1
        {
         sigtot = new float ;    //  total cross-section *JJ
         nhard = new int ;//*JJ number of elementary hard parton-parton scatterings
         npartproj = new int ;//*JJ number of projectile's nucleons participants
         nparttarg = new int ;//*JJ number of target's nucleons participants
         nspecp = new int ;//*JJ number of spectators protons
         nspecn = new int ;//*JJ number of spectators neutrons
        }
        phi = new float ; 
        phir = new float ;
        psi2 = new float ; 
        psi3 = new float ; 
        psi4 = new float ; 
        psi5 = new float ; 
        ecci2 = new float ; 
        ecci3 = new float ; 
        ecci4 = new float ; 
        ecci5 = new float ;  //-----------------

        sprintf(name,"teposevent%i",npom);
        tree[npom] = new TTree(name,"particles"); 
        tree[npom]->Branch("np",np,"np/I");
        tree[npom]->Branch("bim",bim,"bim/F");
        tree[npom]->Branch("nev",nev,"nev/I"); //----new---->
        tree[npom]->Branch("npt",npt,"npt/I");
        tree[npom]->Branch("ngl",ngl,"ngl/I");
        tree[npom]->Branch("kol",kol,"kol/I");
        if(iextree>=1) //extension 1
        {
         tree[npom]->Branch("sigtot",sigtot,"sigtot/F");//*JJ
         tree[npom]->Branch("nhard",nhard,"nhard/I");//*JJ
         tree[npom]->Branch("npartproj",npartproj,"npartproj/I");//*JJ
         tree[npom]->Branch("nparttarg",nparttarg,"nparttarg/I");//*JJ
         tree[npom]->Branch("nspecp",nspecp,"nspecp/I");//*JJ
         tree[npom]->Branch("nspecn",nspecn,"nspecn/I");//*JJ
        }
        tree[npom]->Branch("phi",phi,"phi/F");
        tree[npom]->Branch("phir",phir,"phir/F");
        tree[npom]->Branch("psi2",psi2,"psi2/F");
        tree[npom]->Branch("psi3",psi3,"psi3/F");
        tree[npom]->Branch("psi4",psi4,"psi4/F");
        tree[npom]->Branch("psi5",psi5,"psi5/F");
        tree[npom]->Branch("ecci2",ecci2,"ecci2/F");
        tree[npom]->Branch("ecci3",ecci3,"ecci3/F");
        tree[npom]->Branch("ecci4",ecci4,"ecci4/F");
        tree[npom]->Branch("ecci5",ecci5,"ecci5/F"); //----------- 
        tree[npom]->Branch("zus",zus,"zus[np]/F"); 
        tree[npom]->Branch("px",px,"px[np]/F"); 
        tree[npom]->Branch("py",py,"py[np]/F"); 
        tree[npom]->Branch("pz",pz,"pz[np]/F"); 
        tree[npom]->Branch("e",e,"e[np]/F"); 
        tree[npom]->Branch("x",x,"x[np]/F"); 
        tree[npom]->Branch("y",y,"y[np]/F"); 
        tree[npom]->Branch("z",z,"z[np]/F"); 
        tree[npom]->Branch("t",t,"t[np]/F"); 
        tree[npom]->Branch("id",  id, "id[np]/I"); 
        tree[npom]->Branch("ist",ist,"ist[np]/I"); 
        tree[npom]->Branch("ity",ity,"ity[np]/I"); 
        tree[npom]->Branch("ior",ior,"ior[np]/I"); 
        tree[npom]->Branch("jor",jor,"jor[np]/I"); 
}

extern "C" void treehead_(int *_npom, int *_iversn, int *_laproj, int *_maproj, int *_latarg, int *_matarg, float *_engy, int *_nfull, int *_nfreeze){
        int npom=*_npom-1;
        *iversn = *_iversn ;
        *laproj = *_laproj ;
        *maproj = *_maproj ;
        *latarg = *_latarg ;
        *matarg = *_matarg ;
        *engy   = *_engy ;
        *nfull  = *_nfull ;
        *nfreeze= *_nfreeze ;
        thead[npom]->Fill() ;   
}
 
extern "C" void filltree_(int *_npom, int *_np, float *_bim, float *_sigtot, int *_iextree
           , int *_nev, int *_npt, int *_ngl, int *_kol, int *_nhard, int *_npartproj, int *_nparttarg, int *_nspecp, int *_nspecn //*JJ
           , float *_phi, float *_phir, float *_psi2, float *_psi3, float *_psi4, float *_psi5, float *_ecci2, float *_ecci3, float *_ecci4, float *_ecci5
           , int *_id, int *_ist, int *_ity, int *_ior, int *_jor, float *_zus, float *_px, float *_py, float *_pz,  float *_e
           , float *_x, float *_y, float *_z, float *_t, int *_iret) {
        int npom=*_npom-1;
        int iextree=*_iextree;
        *_iret = 0;
        *np = *_np ;
        if(*np>ntotmax){
          *_iret=1;
          return;
        }
        //if(npom == 0){ std::cout<<" *****fil****  "<<  npom <<"  NP "<<  *np <<std::endl; }
        *bim = *_bim ;
        *nev = *_nev ; //----new---->
        *npt = *_npt ;
        *ngl = *_ngl ;
        *kol = *_kol ;
        if(iextree>=1) //extension 1
        {
        *sigtot = *_sigtot ;//*JJ
        *nhard = *_nhard ;//*JJ
        *npartproj = *_npartproj ;//*JJ
        *nparttarg = *_nparttarg ;//*JJ
        *nspecp = *_nspecp ;//*JJ
        *nspecn = *_nspecn ;//*JJ
        }
        *phi = *_phi ;
        *phir = *_phir ;
        *psi2 = *_psi2 ;
        *psi3 = *_psi3 ;
        *psi4 = *_psi4 ;
        *psi5 = *_psi5 ;
        *ecci2 = *_ecci2 ;
        *ecci3 = *_ecci3 ;
        *ecci4 = *_ecci4 ;
        *ecci5 = *_ecci5 ; //----------- 
        for(int i=0; i<*np; i++){
        id[i] = _id[i] ;
        ist[i] = _ist[i] ;
        ity[i] = _ity[i] ;
        ior[i] = _ior[i]-1 ;
        jor[i] = _jor[i]-1 ;
        zus[i] = _zus[i] ;
        px[i] = _px[i] ;
        py[i] = _py[i] ;
        pz[i] = _pz[i] ;
        e[i] = _e[i] ;
        x[i] = _x[i] ;
        y[i] = _y[i] ;
        z[i] = _z[i] ;
        t[i] = _t[i] ;
        }
        // std::cout << " fill tree with " << *np << " particles" << std::endl ;
        //if(npom == 0){ std::cout<<" *****fil****   Fill " <<std::endl; }
        tree[npom]->Fill() ;    
}

extern "C" void closetree_(int *_npom, char* filename){
//*JJ    std::cout << "******CLOSING THE TREE*******" << std::endl;
       // std::cout << " close tree " << std::endl ;
       int npom=*_npom-1;

       file[npom]->cd() ;       
       //thead[npom]->Write() ;
       //tree[npom]->Write() ;
       file[npom]->Write() ;   //writes all objects into file
       file[npom]->Close() ;
       
       char name [200] ; 
       strcpy(name,filename);
       //std::cout << " closetree ************* " << name << " *************" <<std::endl ;
       TFile *f = new TFile(name,"UPDATE");
       
       sprintf(name,"teposhead%i",npom);
       TTree *th; 
       f->GetObject(name,th); 
       th->SetName("teposhead");
       
       sprintf(name,"teposevent%i",npom);
       TTree *te; 
       f->GetObject(name,te); 
       te->SetName("teposevent");

       f->Write() ;  
       f->Close() ;        
}

