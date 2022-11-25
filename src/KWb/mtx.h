//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

//==============================================================================
//         
//                  Mtx class (Matrix of histograms)
//                  a few centrality class
//
//==============================================================================
#ifndef MTX_H
#define MTX_H
//includes...
#include <iostream>
#include <fstream>
using namespace std;


//==============================================================================
   class mtx
//==============================================================================
{
public:
  //to do:
  //TH1D **hopt;

  //S-func
  TH1D ***hco;
  TH1D ***hcomo;
  TH1D ***hradi;

  TH1D ****hRoslall;
  TH1D ****hRosldir;  

  //CF's
  TH1D ****hcfA;
  TH1D ****hcfo;
  TH1D ****hcfs;
  TH1D ****hcfl;
  TH3D ****h3cf;
  //single technical histograms...
  TH1F **hk_t;
  TH2F **hkt_q;
  TH1F **hNevents;
  TH1F **hmixnum;
  //.......................................................................
  mtx(
      Int_t iCase,             //current case
      Int_t kmin, Int_t kmax,  //min and max centralities are using in job
      Int_t Ce,                //total number of centralities (>=kmax-kmin+1)
      Int_t nkts, Float_t r,
      Int_t bin1, Double_t xm1, 
      Int_t bin3, Double_t xm3) 
  //'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''   
  {
    iCase--;//cases in fortran are 1,2,3,4,... but 0,1,2,3,.. in c++ 
    cout<<"MTX class: iCase "<<iCase<<" NCentralities "<<Ce<<" kmin "<<kmin<<" kmax "<<kmax<<" NkTs "<<nkts<<endl;
    pCentrality=Ce;
    pnCase=iCase;
    pkt=nkts;
    TString cc [10] = {"Kt0","Kt1","Kt2","Kt3","Kt4","Kt5","Kt6","Kt7","Kt8","Kt9"} ;
    TString scase [10] = {"Ca1","Ca2","Ca3","Ca4","Ca5","Ca6","Ca7","Ca8","Ca9","CaA"} ;
    TString xyzt [4] = {"x","y","z","t"};
    TString scentrality[10] = {"Ce1","Ce2","Ce3","Ce4","Ce5","Ce6","Ce7","Ce8","Ce9","CeA"};
    TString scf[3] = {"Qr","Qn","Qm"};//real,norm,mixed
    if(kmax<kmin) {cout<<" Error:  kmax<kmin, kmin="<<kmin<<" kmax="<<kmax<<endl;exit(-1);}
    if(kmax>10)   {cout<<" Error:  kmax>10, kmax="<<kmax<<endl;exit(-1);}
    Int_t uCe=kmax-kmin+1;//number of really used centrality classes
    puCentrality=uCe;
    
    //histograms integrated over kt.....
    //and different centrality class
    hco = new TH1D** [uCe];
    hradi = new TH1D** [uCe];
    hcomo = new TH1D** [uCe];
    hk_t = new TH1F* [uCe];
    hkt_q = new TH2F* [uCe];
    hNevents = new TH1F* [uCe];
    hmixnum = new TH1F* [uCe];
    for(Int_t uc=0; uc<uCe; uc++) {
     Int_t uc1=kmin+uc-1;     //   ****** here -1 necessary !!!!
     hco[uc] = new TH1D* [4] ;
     hradi[uc] = new TH1D* [4] ;
     hcomo[uc] = new TH1D* [4] ;
     TString n1 = "hk_t"+scase[iCase]+scentrality[uc1];
      hk_t[uc] = new TH1F(n1,"dN/dk_{T}"+scase[iCase]+scentrality[uc1],100,0.,1.);
      //
     TString n2 = "hkt_q"+scase[iCase]+scentrality[uc1];
     hkt_q[uc] = new TH2F(n2,"k_{T} vs Qinv"+scase[iCase]+scentrality[uc1],100,0.,1,100,0.,1.);
      //
     TString nNevents = "hNevents"+scase[iCase]+scentrality[uc1];
      hNevents[uc] = new TH1F(nNevents,"Number of selected events"+scase[iCase]+scentrality[uc1],1,1.,2.);
     TString nmixnum = "hmixnum"+scase[iCase]+scentrality[uc1];
      hmixnum[uc] = new TH1F(nmixnum,"Number of part in mixing buffer"+scase[iCase]+scentrality[uc1],1,1.,2.);
     for(Int_t j=0;j<4;j++) {
     TString nco="hR"+xyzt[j]+scase[iCase]+scentrality[uc1];
      hco[uc][j] = new TH1D(nco,nco, 100,-r,r) ;
     TString ncomo="hRp"+xyzt[j]+scase[iCase]+scentrality[uc1];
      hcomo[uc][j] = new TH1D(ncomo,ncomo, 100,-r,r) ;
     TString nradi="hradi"+xyzt[j]+scase[iCase]+scentrality[uc1];
      hradi[uc][j] = new TH1D(nradi,nradi, 100,-r,r) ;
     }
    }//for(Int_t uc=0; uc<uCe; uc++){


    //histograms with kT dependence...
    hcfA = new TH1D*** [uCe];
    hcfo =     new TH1D*** [uCe];
    hcfs =     new TH1D*** [uCe];
    hcfl =     new TH1D*** [uCe];
    h3cf = new TH3D*** [uCe];
    hRoslall = new TH1D*** [uCe];
    hRosldir = new TH1D*** [uCe];
    
    for(Int_t uc=0; uc<uCe; uc++) {
      Int_t uc1=kmin+uc-1;  //   ****** here -1 necessary !!!!;
      hcfA[uc] = new TH1D** [nkts];
      hcfo[uc] = new TH1D** [nkts];
      hcfs[uc] = new TH1D** [nkts];
      hcfl[uc] = new TH1D** [nkts];
      h3cf[uc] = new TH3D** [nkts];

      hRoslall[uc] = new TH1D** [nkts];
      hRosldir[uc] = new TH1D** [nkts];

      for(int j=0; j<nkts; j++){//j=0 is to over all kT
      hcfA[uc][j] = new TH1D* [3];
      hcfo[uc][j] = new TH1D* [3];
      hcfs[uc][j] = new TH1D* [3];
      hcfl[uc][j] = new TH1D* [3];
      h3cf[uc][j] = new TH3D* [3];

      for(Int_t K=0;K<3;K++){
	TString scfA="h"+scf[K]+"a"+scase[iCase]+scentrality[uc1]+cc[j];
	 hcfA[uc][j][K]=new TH1D(scfA,scfA, bin1,0.,xm1);
 	TString scfo="h"+scf[K]+"o"+scase[iCase]+scentrality[uc1]+cc[j];
	 hcfo[uc][j][K]=new TH1D(scfo,scfo, bin1,0.,xm1);
 	TString scfs="h"+scf[K]+"s"+scase[iCase]+scentrality[uc1]+cc[j];
	 hcfs[uc][j][K]=new TH1D(scfs,scfs, bin1,0.,xm1);
 	TString scfl="h"+scf[K]+"l"+scase[iCase]+scentrality[uc1]+cc[j];
	 hcfl[uc][j][K]=new TH1D(scfl,scfl, bin1,0.,xm1);
 	TString s3cf="h3"+scf[K]+scase[iCase]+scentrality[uc1]+cc[j];
	 h3cf[uc][j][K]=new TH3D(s3cf,s3cf, bin3,0.,xm3,bin3,0.,xm3,bin3,0.,xm3);
     }

      hRoslall[uc][j]= new TH1D* [4] ;
      hRosldir[uc][j]= new TH1D* [4] ;
      for(int m=0;m<4;m++) {
	TString nRoslall="hRall"+xyzt[m]+scase[iCase]+scentrality[uc1]+cc[j];
	hRoslall[uc][j][m] = new TH1D(nRoslall,nRoslall, 100,-r,r) ;
	TString nRosldir="hRdir"+xyzt[m]+scase[iCase]+scentrality[uc1]+cc[j];
	hRosldir[uc][j][m] = new TH1D(nRosldir,nRosldir, 100,-r,r) ;
      }
      }     
    }
    
    
  };
  //...............................................................
  void WriteHistos(TString OutHBTHisto,Int_t iCase, Bool_t bUpdate)
  //'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  {
    TFile *of;
    //bUpdate=false;//write at the end of all cases
    //cout<<OutHBTHisto<<" Bool_t bUpdate= "<<bUpdate<<endl;
    if(bUpdate) {
      of = new TFile(OutHBTHisto, "UPDATE");
      cout<<"Write matrix histograms to updated file: "<<OutHBTHisto<<endl;
    }
    else { 
      of = new TFile(OutHBTHisto, "RECREATE");
      cout<<"Write matrix histograms to the new file: "<<OutHBTHisto<<endl;
    }
    //cout <<"Write Matrix histogram "<<pkt<<endl;
    // Write out result histograms
    for(Int_t uc=0; uc<puCentrality; uc++) {
      //single technical histograms...
      hk_t[uc]->Write();
      hkt_q[uc]->Write();
      hNevents[uc]->Write();
      hmixnum[uc]->Write();
      for(Int_t j=0;j<4;j++){
	hco[uc][j]->Write();
	hcomo[uc][j]->Write();
	hradi[uc][j]->Write();
      }
    }
    for(Int_t uc=0; uc<puCentrality; uc++) {
      for(Int_t j=0;j<pkt;j++){
	for(Int_t K=0;K<3;K++){
	hcfA[uc][j][K]->Write();
	hcfo[uc][j][K]->Write();
	hcfs[uc][j][K]->Write();
	hcfl[uc][j][K]->Write();
	h3cf[uc][j][K]->Write();
	}
	for(int m=0;m<4;m++) {
	  hRoslall[uc][j][m]->Write();
	  hRosldir[uc][j][m]->Write();
	}
      }
    }

    of->Close();
  };
  //..............
  virtual ~mtx(){
  //''''''''''''''
    cout <<"MTX class: Delete Matrix histogram "<<endl;
    for(Int_t uc=0; uc<puCentrality; uc++) {
      //single technical histograms...
      hk_t[uc]->Delete();
      hNevents[uc]->Delete();
      hmixnum[uc]->Delete();
      for(Int_t j=0;j<4;j++){
	hco[uc][j]->Delete();
	hcomo[uc][j]->Delete();
	hradi[uc][j]->Delete();
      }
    }

    for(Int_t uc=0; uc<puCentrality; uc++) {
      for(Int_t j=0;j<pkt;j++){
	for(Int_t K=0;K<3;K++){
	hcfA[uc][j][K]->Delete();
	hcfo[uc][j][K]->Delete();
	hcfs[uc][j][K]->Delete();
	hcfl[uc][j][K]->Delete();
	h3cf[uc][j][K]->Delete();
	}
	for(Int_t m=0;m<4;m++) {
	  hRoslall[uc][j][m]->Delete();
	  hRosldir[uc][j][m]->Delete();
	}
      }
    }

  };

protected:
  Int_t pCentrality;
  Int_t puCentrality;
  Int_t pnCase;
  Int_t pkt;

#ifdef __ROOT__
  ClassDef(EpMaHi,1)
#endif

};

#endif
