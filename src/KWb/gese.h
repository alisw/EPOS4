//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

//----------------------------------------------------------//
//                                                          //
// Class to do get/set  epos and                            //
//                      fsi parameters for Lednicky's code  //
//                                                          //
//                                                          //
// author: Konstantin Mikhailov (kmikhail@itep.ru)          //
//----------------------------------------------------------//

#ifndef GETTERSETTER_H
#define GETTERSETTER_H
//cout, endl,...
#include <iostream>
#include <fstream>
using namespace std;

class gettersetter
{

public:
  gettersetter(){};
  //Getters
  //opt...
  Int_t getnout            ()  const { return nout_              ;}
  Int_t getnClass          ()  const { return nClass_            ;}
  Float_t getbimpMin       ()  const { return bimpMin_           ;}
  Float_t getbimpMax       ()  const { return bimpMax_           ;}
  Int_t getnKtBins         ()  const { return nKtBins_           ;}
  Float_t getKtMin         ()  const { return KtMin_             ;}
  Float_t getKtMax         ()  const { return KtMax_             ;}
  Int_t getnCase           ()  const { return nCases_             ;}
  Int_t getiCase           ()  const { return iCase_             ;}
  Int_t getMixedEvents     ()  const { return MixedEvents_       ;}
  Int_t getMaxNpart        ()  const { return MaxNpart_          ;}
  //fsi...
  Int_t getPa1            ()  const { return Pa1_              ;}
  Int_t getPa2            ()  const { return Pa2_              ;} 
  Int_t getLL             ()  const { return LL_               ;}
  Int_t getITEST          ()  const { return ITEST_            ;}
  Int_t getICH            ()  const { return  ICH_             ;}        
  Int_t getIQS            ()  const { return  IQS_             ;}        
  Int_t getISI            ()  const { return  ISI_	       ;}
  Int_t getI3C            ()  const { return  I3C_	       ;}
  Int_t getNUCLMASS       ()  const { return  NUCLMASS_        ;}
  Int_t getNUCLCHARGE     ()  const { return  NUCLCHARGE_      ;}
  Int_t getNUCLCHARGESIGN ()  const { return  NUCLCHARGESIGN_  ;}
  Int_t getSphereApp      ()  const { return  SphereApp_       ;}	     
  Int_t getT0App          ()  const { return  T0App_	       ;}
  Int_t getfastQS         ()  const { return  fastQS_          ;} 
  //real iOpt(10)
  Float_t getiOpt0       () const { return iOpt0_             ;} 
  Float_t getiOpt1       () const { return iOpt1_             ;} 
  Float_t getiOpt2       () const { return iOpt2_             ;} 
  Float_t getiOpt3       () const { return iOpt3_             ;} 
  Float_t getiOpt4       () const { return iOpt4_             ;} 
  Float_t getiOpt5       () const { return iOpt5_             ;} 
  Float_t getiOpt6       () const { return iOpt6_             ;} 
  Float_t getiOpt7       () const { return iOpt7_             ;} 
  Float_t getiOpt8       () const { return iOpt8_             ;} 
  Float_t getiOpt9       () const { return iOpt9_             ;} 
  //real dcut(20)
  Float_t getdcut00() const { return dcut00_;}
  Float_t getdcut01() const { return dcut01_;}
  Float_t getdcut02() const { return dcut02_;}
  Float_t getdcut03() const { return dcut03_;}
  Float_t getdcut04() const { return dcut04_;}
  Float_t getdcut05() const { return dcut05_;}
  Float_t getdcut06() const { return dcut06_;}
  Float_t getdcut07() const { return dcut07_;}
  Float_t getdcut08() const { return dcut08_;}
  Float_t getdcut09() const { return dcut09_;}
  Float_t getdcut10() const { return dcut10_;}
  Float_t getdcut11() const { return dcut11_;}
  Float_t getdcut12() const { return dcut12_;}
  Float_t getdcut13() const { return dcut13_;}
  Float_t getdcut14() const { return dcut14_;}
  Float_t getdcut15() const { return dcut15_;}
  Float_t getdcut16() const { return dcut16_;}
  Float_t getdcut17() const { return dcut17_;}
  Float_t getdcut18() const { return dcut18_;}
  Float_t getdcut19() const { return dcut19_;}
  //char file names etc...
  //char* getfnamein          ()   const { return fnamein_           ;}
  char* getfnamein          ()   { return fnamein_           ;}
  //char file names etc...
  char* getcbasout          ()   { return cbasout_           ;}



  //Setters 
  //real iOpt(10)
  void setiOpt0     (Float_t i0) { iOpt0_ = i0; }
  void setiOpt1     (Float_t i1) { iOpt1_ = i1; }
  void setiOpt2     (Float_t i2) { iOpt2_ = i2; }
  void setiOpt3     (Float_t i3) { iOpt3_ = i3; }
  void setiOpt4     (Float_t i4) { iOpt4_ = i4; }
  void setiOpt5     (Float_t i5) { iOpt5_ = i5; }
  void setiOpt6     (Float_t i6) { iOpt6_ = i6; }
  void setiOpt7     (Float_t i7) { iOpt7_ = i7; }
  void setiOpt8     (Float_t i8) { iOpt8_ = i8; }
  void setiOpt9     (Float_t i9) { iOpt9_ = i9; }
  //real dcut(20)
  void setdcut00  (Float_t d00) {dcut00_ = d00;}
  void setdcut01  (Float_t d01) {dcut01_ = d01;}
  void setdcut02  (Float_t d02) {dcut02_ = d02;}
  void setdcut03  (Float_t d03) {dcut03_ = d03;}
  void setdcut04  (Float_t d04) {dcut04_ = d04;}
  void setdcut05  (Float_t d05) {dcut05_ = d05;}
  void setdcut06  (Float_t d06) {dcut06_ = d06;}
  void setdcut07  (Float_t d07) {dcut07_ = d07;}
  void setdcut08  (Float_t d08) {dcut08_ = d08;}
  void setdcut09  (Float_t d09) {dcut09_ = d09;}
  void setdcut10  (Float_t d10) {dcut10_ = d10;}
  void setdcut11  (Float_t d11) {dcut11_ = d11;}
  void setdcut12  (Float_t d12) {dcut12_ = d12;}
  void setdcut13  (Float_t d13) {dcut13_ = d13;}
  void setdcut14  (Float_t d14) {dcut14_ = d14;}
  void setdcut15  (Float_t d15) {dcut15_ = d15;}
  void setdcut16  (Float_t d16) {dcut16_ = d16;}
  void setdcut17  (Float_t d17) {dcut17_ = d17;}
  void setdcut18  (Float_t d18) {dcut18_ = d18;}
  void setdcut19  (Float_t d19) {dcut19_ = d19;}

  //char file names etc...
  void setfnamein         (char *ca1)  { 
    sprintf(fnamein_,"%s",ca1);
  }
  void setcbasout         (char *ca2)  { 
    sprintf(cbasout_,"%s",ca2);
  }
  //opt...
  void setnout            (Int_t o0)   { nout_               = o0;}
  void setnClass          (Int_t o1)   { nClass_             = o1;}
  void setbimpMin         (Float_t o2) { bimpMin_            = o2;}
  void setbimpMax         (Float_t o3) { bimpMax_            = o3;}
  void setnKtBins         (Int_t o4)   { nKtBins_            = o4;}
  void setKtMin           (Float_t o5) { KtMin_              = o5;}
  void setKtMax           (Float_t o6) { KtMax_              = o6;}
  void setnCases           (Int_t o71) { nCases_            = o71;}
  void setiCase           (Int_t o72)  { iCase_             = o72;}
  void setMixedEvents     (Int_t o7)   { MixedEvents_        = o7;}
  void setMaxNpart        (Int_t o8)   { MaxNpart_           = o8;}
  //fsi...
  void setPa1            (Int_t h1) { Pa1_            = h1; }
  void setPa2            (Int_t h2) { Pa2_            = h2; }
  void setLL             (Int_t h3) { LL_             = h3; }
  void setITEST          (Int_t h4) { ITEST_          = h4; }
  void setICH            (Int_t h5) { ICH_            = h5; }          
  void setIQS            (Int_t h6) { IQS_            = h6; }          
  void setISI	         (Int_t h7) { ISI_	      = h7; }
  void setI3C	         (Int_t h8) { I3C_	      = h8; }
  void setNUCLMASS       (Int_t h9) { NUCLMASS_       = h9; } 
  void setNUCLCHARGE     (Int_t hq) { NUCLCHARGE_     = hq; } 
  void setNUCLCHARGESIGN (Int_t hw) { NUCLCHARGESIGN_ = hw; } 
  void setSphereApp      (Int_t he) { SphereApp_      = he; }	     
  void setT0App	         (Int_t hr) { T0App_	      = hr; }    
  void setfastQS         (Int_t ht) { fastQS_         = ht; }   
   

  virtual   ~gettersetter(){};

private:

  //char file names etc...
  char    fnamein_[150]  ;
  char    cbasout_[150]  ;
  //opt...
  Int_t    nout_         ;  
  Int_t    nClass_       ;  
  Float_t  bimpMin_      ; 
  Float_t  bimpMax_      ; 
  Int_t    nKtBins_      ; 
  Float_t  KtMin_        ; 
  Float_t  KtMax_        ; 
  Int_t    nCases_  ; 
  Int_t    iCase_  ; 
  Int_t    MixedEvents_  ; 
  Int_t	   MaxNpart_     ;
  //real iOpt(10)
  Float_t iOpt0_;
  Float_t iOpt1_;
  Float_t iOpt2_;
  Float_t iOpt3_;
  Float_t iOpt4_;
  Float_t iOpt5_;
  Float_t iOpt6_;
  Float_t iOpt7_;
  Float_t iOpt8_;
  Float_t iOpt9_;
  //real dcut(20)
  Float_t dcut00_;
  Float_t dcut01_;
  Float_t dcut02_;
  Float_t dcut03_;
  Float_t dcut04_;
  Float_t dcut05_;
  Float_t dcut06_;
  Float_t dcut07_;
  Float_t dcut08_;
  Float_t dcut09_;
  Float_t dcut10_;
  Float_t dcut11_;
  Float_t dcut12_;
  Float_t dcut13_;
  Float_t dcut14_;
  Float_t dcut15_;
  Float_t dcut16_;
  Float_t dcut17_;
  Float_t dcut18_;
  Float_t dcut19_;


  //fsi...
  Int_t Pa1_;
  Int_t Pa2_;
  Int_t LL_;
  Int_t ITEST_;
  Int_t ICH_;            
  Int_t IQS_;            
  Int_t ISI_;	     
  Int_t I3C_;	     
  Int_t NUCLMASS_;       
  Int_t NUCLCHARGE_;     
  Int_t NUCLCHARGESIGN_; 
  Int_t SphereApp_;	     
  Int_t T0App_;	     
  Int_t fastQS_;       
  
#ifdef __ROOT__
  ClassDef(gettersetter,1)
#endif
};


#endif
