//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>

#include "ihi.h"

TH2D **h1, **h2;
TFile *outputFile ;

int kkmax ; 

using namespace std ;

void open2dhisto_(char * filename, int * kmax, int * nxbins, int * nybins, double * etamax)
{
	outputFile = new TFile(filename, "RECREATE");
	outputFile->cd();
	
	kkmax = *kmax ;
	
	h1 = new TH2D* [kkmax] ;
	char name1 [10][10] = {"1K1","1K2","1K3","1K4","1K5","1K6","1K7","1K8","1K9","1K10"} ;
	char title1 [10][10] = {"1K1","1K2","1K3","1K4","1K5","1K6","1K7","1K8","1K9","1K10"} ;
	for(int i=0; i<kkmax; i++){
	h1[i] = new TH2D(name1[i], title1[i], *nxbins, -*etamax, *etamax, *nybins, -1.0*TMath::Pi(), TMath::Pi() ) ;
	}

	h2 = new TH2D* [kkmax] ;
	char name2 [10][10] = {"2K1","2K2","2K3","2K4","2K5","2K6","2K7","2K8","2K9","2K10"} ;
	char title2 [10][10] = {"2K1","2K2","2K3","2K4","2K5","2K6","2K7","2K8","2K9","2K10"} ;
	for(int i=0; i<kkmax; i++){
	h2[i] = new TH2D(name2[i], title2[i], *nxbins, -*etamax, *etamax, *nybins, -1.0*TMath::Pi(), TMath::Pi() ) ;
	}
	
}

void fill2dhisto1_(int * kk, double * delta_eta, double * delta_phi)
{
	//cout << *delta_eta << " " << *delta_phi << endl;
        h1[*kk-1]->Fill(*delta_eta, *delta_phi) ;
}

void fill2dhisto2_(int * kk, double * delta_eta, double * delta_phi)
{
	h2[*kk-1]->Fill(*delta_eta, *delta_phi) ;
}



void close2dhisto_()
{
	// call writing methods for histograms instead of writing everything into current ROOT file
	for(int i=0; i<kkmax; i++){
        outputFile->cd();
	h1[i]->Write() ;
	h2[i]->Write() ;
	}
//	outputFile->Write() ;
	outputFile->Close() ;
}
