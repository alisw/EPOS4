//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

class ETree{
public:
  ETree(){};
  TTree *EposHeadTree;          
  TTree *EposTree;       
  TFile* fEposTreeFile;
  ETree(char* fnamein){
    //cout<<"---------------- Read Epos Tree file ------------------"<<endl;
    cout<<"ETree: Open " << fnamein <<endl;
    fEposTreeFile = TFile::Open(fnamein); 
    //cout<<"  file "<<fEposTreeFile<<endl;
    if(fEposTreeFile==NULL) {
      cout<<"Error: can't open file "<< fnamein <<endl;
      exit(-1);
    } 
    TTree * hd=(TTree*)fEposTreeFile->Get("teposhead");
    EposHeadTree=(TTree*)hd;
    TTree * ep=(TTree*)fEposTreeFile->Get("teposevent");
    EposTree=(TTree*)ep;           
    //cout<<"------------- End of Read Epos Tree file --------------"<<endl;
  };
  ~ETree(){fEposTreeFile->TFile::Close(); };
private:
protected:
};


