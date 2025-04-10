//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "ak.h"

using namespace std ;

//-------------------------------------------------------------------------------------------

Mudiar<float> *dmass_c0 ;
extern "C" void dmass_c0create_(int* n1)         { dmass_c0 = new Mudiar<float>(*n1) ;  }
extern "C" void dmass_c0destroy_(void)           { delete dmass_c0 ; }
extern "C" void dmass_c0set_(int* i1,float* val) { dmass_c0->set(*i1-1, *val) ;  }
extern "C" void dmass_c0get_(int* i1,float* val) { dmass_c0->get(*i1-1, *val) ;  }
extern "C" void dmass_c0acc_(int* i1,float* val) { dmass_c0->acc(*i1-1, *val) ;  }
Mudiar<float> *dmass_c1 ;
extern "C" void dmass_c1create_(int* n1,int* n2)         { dmass_c1 = new Mudiar<float>(*n1,*n2) ;  }
extern "C" void dmass_c1destroy_(void)                   { delete dmass_c1 ; }
extern "C" void dmass_c1set_(int* i1,int* i2,float* val) { dmass_c1->set(*i1-1,*i2-1, *val) ;  }
extern "C" void dmass_c1get_(int* i1,int* i2,float* val) { dmass_c1->get(*i1-1,*i2-1, *val) ;  }
extern "C" void dmass_c1acc_(int* i1,int* i2,float* val) { dmass_c1->acc(*i1-1,*i2-1, *val) ;  }
Mudiar<float> *dmass_c2 ;
extern "C" void dmass_c2create_(int* n1,int* n2, int* n3)          { dmass_c2 = new Mudiar<float>(*n1,*n2,*n3) ;  }
extern "C" void dmass_c2destroy_(void)                             { delete dmass_c2 ; }
extern "C" void dmass_c2set_(int* i1, int* i2, int* i3,float* val) { dmass_c2->set(*i1-1, *i2-1, *i3-1, *val) ;  }
extern "C" void dmass_c2get_(int* i1, int* i2, int* i3,float* val) { dmass_c2->get(*i1-1, *i2-1, *i3-1, *val) ;  }
extern "C" void dmass_c2acc_(int* i1, int* i2, int* i3,float* val) { dmass_c2->acc(*i1-1, *i2-1, *i3-1, *val) ;  }
Mudiar<float> *dmass_c3 ;
extern "C" void dmass_c3create_(int* n1,int* n2, int* n3, int* n4)          { dmass_c3 = new Mudiar<float>(*n1,*n2,*n3,*n4) ;  }
extern "C" void dmass_c3destroy_(void)                                      { delete dmass_c3 ; }
extern "C" void dmass_c3set_(int* i1, int* i2, int* i3, int* i4,float* val) { dmass_c3->set(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void dmass_c3get_(int* i1, int* i2, int* i3, int* i4,float* val) { dmass_c3->get(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void dmass_c3acc_(int* i1, int* i2, int* i3, int* i4,float* val) { dmass_c3->acc(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
Mudiar<float> *dmass_i0 ;
extern "C" void dmass_i0create_(int* n1)         { dmass_i0 = new Mudiar<float>(*n1) ;  }
extern "C" void dmass_i0destroy_(void)           { delete dmass_i0 ; }
extern "C" void dmass_i0set_(int* i1,float* val) { dmass_i0->set(*i1-1, *val) ;  }
extern "C" void dmass_i0get_(int* i1,float* val) { dmass_i0->get(*i1-1, *val) ;  }
extern "C" void dmass_i0acc_(int* i1,float* val) { dmass_i0->acc(*i1-1, *val) ;  }
Mudiar<float> *dmass_i1 ;
extern "C" void dmass_i1create_(int* n1,int* n2)         { dmass_i1 = new Mudiar<float>(*n1,*n2) ;  }
extern "C" void dmass_i1destroy_(void)                   { delete dmass_i1 ; }
extern "C" void dmass_i1set_(int* i1,int* i2,float* val) { dmass_i1->set(*i1-1,*i2-1, *val) ;  }
extern "C" void dmass_i1get_(int* i1,int* i2,float* val) { dmass_i1->get(*i1-1,*i2-1, *val) ;  }
extern "C" void dmass_i1acc_(int* i1,int* i2,float* val) { dmass_i1->acc(*i1-1,*i2-1, *val) ;  }
Mudiar<float> *dmass_i2 ;
extern "C" void dmass_i2create_(int* n1,int* n2, int* n3)          { dmass_i2 = new Mudiar<float>(*n1,*n2,*n3) ;  }
extern "C" void dmass_i2destroy_(void)                             { delete dmass_i2 ; }
extern "C" void dmass_i2set_(int* i1, int* i2, int* i3,float* val) { dmass_i2->set(*i1-1, *i2-1, *i3-1, *val) ;  }
extern "C" void dmass_i2get_(int* i1, int* i2, int* i3,float* val) { dmass_i2->get(*i1-1, *i2-1, *i3-1, *val) ;  }
extern "C" void dmass_i2acc_(int* i1, int* i2, int* i3,float* val) { dmass_i2->acc(*i1-1, *i2-1, *i3-1, *val) ;  }
Mudiar<float> *dmass_i3 ;
extern "C" void dmass_i3create_(int* n1,int* n2, int* n3, int* n4)          { dmass_i3 = new Mudiar<float>(*n1,*n2,*n3,*n4) ;  }
extern "C" void dmass_i3destroy_(void)                                      { delete dmass_i3 ; }
extern "C" void dmass_i3set_(int* i1, int* i2, int* i3, int* i4,float* val) { dmass_i3->set(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void dmass_i3get_(int* i1, int* i2, int* i3, int* i4,float* val) { dmass_i3->get(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void dmass_i3acc_(int* i1, int* i2, int* i3, int* i4,float* val) { dmass_i3->acc(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }

//###################################################################################################
//###################################################################################################

Mudiar<int> *deformopt ;
extern "C" void deformoptcreate_(int* n1)      { deformopt = new Mudiar<int>(*n1) ;  }
extern "C" void deformoptdestroy_(void)        { delete deformopt ; }
extern "C" void deformoptset_(int* i1,int* val){ deformopt->set(*i1-1, *val) ;  }
extern "C" void deformoptget_(int* i1,int* val){ deformopt->get(*i1-1, *val) ;  }


//###################################################################################################
//##################################### bor objects #################################################
//###################################################################################################

Mudiar<float> *xpprbor ;
extern "C" void xpprbor_create_(int* n1,int* n2, int* n3)          { xpprbor = new Mudiar<float>(*n1,*n2,*n3) ;  }
extern "C" void xpprbor_destroy_(void)                             { delete xpprbor ; }
extern "C" void xpprbor_set_(int* i1, int* i2, int* i3,float* val) { xpprbor->set(*i1-1, *i2-1, *i3-1, *val) ;  }
extern "C" void xpprbor_get_(int* i1, int* i2, int* i3,float* val) { xpprbor->get(*i1-1, *i2-1, *i3-1, *val) ;  }

Mudiar<float> *xmprbor ;
extern "C" void xmprbor_create_(int* n1,int* n2, int* n3)          { xmprbor = new Mudiar<float>(*n1,*n2,*n3) ;  }
extern "C" void xmprbor_destroy_(void)                             { delete xmprbor ; }
extern "C" void xmprbor_set_(int* i1, int* i2, int* i3,float* val) { xmprbor->set(*i1-1, *i2-1, *i3-1, *val) ;  }
extern "C" void xmprbor_get_(int* i1, int* i2, int* i3,float* val) { xmprbor->get(*i1-1, *i2-1, *i3-1, *val) ;  }

Mudiar<float> *ptprboo ;
extern "C" void ptprboo_create_(int* n1,int* n2, int* n3, int* n4)          { ptprboo = new Mudiar<float>(*n1,*n2,*n3,*n4) ;  }
extern "C" void ptprboo_destroy_(void)                                      { delete ptprboo ; }
extern "C" void ptprboo_set_(int* i1, int* i2, int* i3, int* i4,float* val) { ptprboo->set(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void ptprboo_get_(int* i1, int* i2, int* i3, int* i4,float* val) { ptprboo->get(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }

Mudiar<float> *rapprboo ;
extern "C" void rapprboo_create_(int* n1,int* n2, int* n3, int* n4)          { rapprboo = new Mudiar<float>(*n1,*n2,*n3,*n4) ;  }
extern "C" void rapprboo_destroy_(void)                                      { delete rapprboo ; }
extern "C" void rapprboo_set_(int* i1, int* i2, int* i3, int* i4,float* val) { rapprboo->set(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void rapprboo_get_(int* i1, int* i2, int* i3, int* i4,float* val) { rapprboo->get(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }

Mudiar<float> *gbpom ;
extern "C" void gbpom_create_(int* n1,int* n2, int* n3, int* n4)          { gbpom = new Mudiar<float>(*n1,*n2,*n3,*n4) ;  }
extern "C" void gbpom_destroy_(void)                                      { delete gbpom ; }
extern "C" void gbpom_set_(int* i1, int* i2, int* i3, int* i4,float* val) { gbpom->set(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void gbpom_get_(int* i1, int* i2, int* i3, int* i4,float* val) { gbpom->get(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }

Mudiar<int> *idbor ;
extern "C" void idbor_create_(int* n1,int* n2, int* n3, int* n4)          { idbor = new Mudiar<int>(*n1,*n2,*n3,*n4) ;  }
extern "C" void idbor_destroy_(void)                                      { delete idbor ; }
extern "C" void idbor_set_(int* i1, int* i2, int* i3, int* i4,int* val) { idbor->set(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }
extern "C" void idbor_get_(int* i1, int* i2, int* i3, int* i4,int* val) { idbor->get(*i1-1, *i2-1, *i3-1, *i4-1, *val) ;  }

//###################################################################################################
//###################################### parameters objects #########################################
//###################################################################################################

//-------------------------------------------------------------------------------------------

Mudiar<float> *foparam ;
extern "C" void foparamcreate_(int* n1)      { foparam = new Mudiar<float>(*n1) ;  }
extern "C" void foparamdestroy_(void)        { delete foparam ; }
extern "C" void foparamset_(int* i1,float* val){ foparam->set(*i1-1, *val) ;  }
extern "C" void foparamget_(int* i1,float* val){ foparam->get(*i1-1, *val) ;  }

//-------------------------------------------------------------------------------------------

Mudiar<float> *remn1param ;
extern "C" void remn1paramcreate_(int* n1)         { remn1param = new Mudiar<float>(*n1) ;  }
extern "C" void remn1paramdestroy_(void)           { delete remn1param ; }
extern "C" void remn1paramset_(int* i1,float* val) { remn1param->set(*i1-1, *val) ;  }
extern "C" void remn1paramget_(int* i1,float* val) { remn1param->get(*i1-1, *val) ;  }
extern "C" void remn1paramset6_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6) { remn1param->set6(*v1,*v2,*v3,*v4,*v5,*v6) ;  }
extern "C" void remn1paramget6_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6) { remn1param->get6(*v1,*v2,*v3,*v4,*v5,*v6) ;  }

//-------------------------------------------------------------------------------------------

Mudiar<float> *pfe3param ;
extern "C" void pfe3paramcreate_(int* n1)         { pfe3param = new Mudiar<float>(*n1) ;  }
extern "C" void pfe3paramdestroy_(void)           { delete pfe3param ; }
extern "C" void pfe3paramset_(int* i1,float* val) { pfe3param->set(*i1-1, *val) ;  }
extern "C" void pfe3paramget_(int* i1,float* val) { pfe3param->get(*i1-1, *val) ;  }
extern "C" void pfe3paramset2_(float* v1,float* v2)  { pfe3param->set2(*v1,*v2) ;  }
extern "C" void pfe3paramget2_(float* v1,float* v2)  { pfe3param->get2(*v1,*v2) ;  }
extern "C" void pfe3paramset3_(float* v1,float* v2,float* v3) { pfe3param->set3(*v1,*v2,*v3) ;  }
extern "C" void pfe3paramget3_(float* v1,float* v2,float* v3) { pfe3param->get3(*v1,*v2,*v3) ;  }
extern "C" void pfe3paramset5_(float* v1,float* v2,float* v3,float* v4,float* v5) { pfe3param->set5(*v1,*v2,*v3,*v4,*v5) ;  }
extern "C" void pfe3paramget5_(float* v1,float* v2,float* v3,float* v4,float* v5) { pfe3param->get5(*v1,*v2,*v3,*v4,*v5) ;  }



Mudiar<float> *pfe6param ;
extern "C" void pfe6paramcreate_(int* n1)         { pfe6param = new Mudiar<float>(*n1) ;  }
extern "C" void pfe6paramdestroy_(void)           { delete pfe6param ; }
extern "C" void pfe6paramset_(int* i1,float* val) { pfe6param->set(*i1-1, *val) ;  }
extern "C" void pfe6paramget_(int* i1,float* val) { pfe6param->get(*i1-1, *val) ;  }
extern "C" void pfe6paramset2_(float* v1,float* v2)  { pfe6param->set2(*v1,*v2) ;  }
extern "C" void pfe6paramget2_(float* v1,float* v2)  { pfe6param->get2(*v1,*v2) ;  }

//-------------------------------------------------------------------------------------------

Mudiar<float> *screen1param ;
extern "C" void screen1paramcreate_(int* n1)         { screen1param = new Mudiar<float>(*n1) ;  }
extern "C" void screen1paramdestroy_(void)           { delete screen1param ; }
extern "C" void screen1paramset_(int* i1,float* val) { screen1param->set(*i1-1, *val) ;  }
extern "C" void screen1paramget_(int* i1,float* val) { screen1param->get(*i1-1, *val) ;  }
extern "C" void screen1paramset6_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6) { screen1param->set6(*v1,*v2,*v3,*v4,*v5,*v6) ;  }
extern "C" void screen1paramget6_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6) { screen1param->get6(*v1,*v2,*v3,*v4,*v5,*v6) ;  }

Mudiar<float> *screen2param ;
extern "C" void screen2paramcreate_(int* n1)         { screen2param = new Mudiar<float>(*n1) ;  }
extern "C" void screen2paramdestroy_(void)           { delete screen2param ; }
extern "C" void screen2paramset_(int* i1,float* val) { screen2param->set(*i1-1, *val) ;  }
extern "C" void screen2paramget_(int* i1,float* val) { screen2param->get(*i1-1, *val) ;  }
extern "C" void screen2paramzero_()                  { screen2param->zero() ;  }
extern "C" void screen2paramset5_(float* v1,float* v2,float* v3,float* v4,float* v5) { screen2param->set5(*v1,*v2,*v3,*v4,*v5) ;  }
extern "C" void screen2paramget5_(float* v1,float* v2,float* v3,float* v4,float* v5) { screen2param->get5(*v1,*v2,*v3,*v4,*v5) ;  }

//-------------------------------------------------------------------------------------------

Mudiar<float> *pdfparam ;
extern "C" void pdfparamcreate_(int* n1)      { pdfparam = new Mudiar<float>(*n1) ;  }
extern "C" void pdfparamdestroy_(void)        { delete pdfparam ; }
extern "C" void pdfparamset_(int* i1,float* val){ pdfparam->set(*i1-1, *val) ;  }
extern "C" void pdfparamget_(int* i1,float* val){ pdfparam->get(*i1-1, *val) ;  }

//-------------------------------------------------------------------------------------------

Mudiar<float> *saturparam ;
extern "C" void saturparamcreate_(int* n1)      { saturparam = new Mudiar<float>(*n1) ;  }
extern "C" void saturparamdestroy_(void)        { delete saturparam ; }
extern "C" void saturparamset_(int* i1,float* val){ saturparam->set(*i1-1, *val) ;  }
extern "C" void saturparamget_(int* i1,float* val){ saturparam->get(*i1-1, *val) ;  }
extern "C" void saturparamset2_(float* v1,float* v2)  { saturparam->set2(*v1,*v2) ;  }
extern "C" void saturparamget2_(float* v1,float* v2)  { saturparam->get2(*v1,*v2) ;  }

//-------------------------------------------------------------------------------------------

Mudiar<float> *core1param ;
extern "C" void core1paramcreate_(int* n1)      { core1param = new Mudiar<float>(*n1) ;  }
extern "C" void core1paramdestroy_(void)        { delete core1param ; }
extern "C" void core1paramset_(int* i1,float* val){ core1param->set(*i1-1, *val) ;  }
extern "C" void core1paramget_(int* i1,float* val){ core1param->get(*i1-1, *val) ;  }
extern "C" void core1paramset4_(float* v1,float* v2,float* v3,float* v4) { core1param->set4(*v1,*v2,*v3,*v4) ;  }
extern "C" void core1paramget4_(float* v1,float* v2,float* v3,float* v4) { core1param->get4(*v1,*v2,*v3,*v4) ;  }
extern "C" void core1paramset7_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6,float* v7) { core1param->set7(*v1,*v2,*v3,*v4,*v5,*v6,*v7) ;  }
extern "C" void core1paramget7_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6,float* v7) { core1param->get7(*v1,*v2,*v3,*v4,*v5,*v6,*v7) ;  }
Mudiar<float> *feloss ;
extern "C" void felosscreate_(int* n1)      { feloss = new Mudiar<float>(*n1) ;  }
extern "C" void felossdestroy_(void)        { delete feloss ; }
extern "C" void felossset_(int* i1,float* val){ feloss->set(*i1-1, *val) ;  }
extern "C" void felossget_(int* i1,float* val){ feloss->get(*i1-1, *val) ;  }
extern "C" void felossset6_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6) { feloss->set6(*v1,*v2,*v3,*v4,*v5,*v6) ;  }
extern "C" void felossget6_(float* v1,float* v2,float* v3,float* v4,float* v5,float* v6) { feloss->get6(*v1,*v2,*v3,*v4,*v5,*v6) ;  }

//###################################################################################################
//################################ maxsize object -> dimensions #####################################
//###################################################################################################

Mudiar<int> *maxsize ;
extern "C" void maxsize_create_(int* n1)      { maxsize = new Mudiar<int>(*n1) ;  }
extern "C" void maxsize_destroy_(void)        { delete maxsize ; }
extern "C" void maxsize_set_(int* i1,int* val){ maxsize->set(*i1-1, *val) ;  }
extern "C" void maxsize_get_(int* i1,int* val){ maxsize->get(*i1-1, *val) ;  }

//###################################################################################################
//################################### event variables object ########################################
//###################################################################################################

//-------------------------------------------------------------------------------------------

Mudiar<float> *eventvari ;
extern "C" void eventvaricreate_(int* n1)         { eventvari = new Mudiar<float>(*n1) ;  }
extern "C" void eventvaridestroy_(void)           { delete eventvari ; }
extern "C" void eventvariset_(int* i1,float* val) { eventvari->set(*i1-1, *val) ;  }
extern "C" void eventvariget_(int* i1,float* val) { eventvari->get(*i1-1, *val) ;  }

//###################################################################################################
//###################################################################################################

//-------------------------------------------------------------------------------------------

Mudiar<int> *lspecs ;
extern "C" void lspecscreate_(int* n1,int* n2)      { lspecs = new Mudiar<int>(*n1,*n2) ;  }
extern "C" void lspecsdestroy_(void)                { delete lspecs ; }
extern "C" void lspecsset_(int* i1,int* i2,int* val){ lspecs->set(*i1-1,*i2-1, *val) ;  }
extern "C" void lspecsget_(int* i1,int* i2,int* val){ lspecs->get(*i1-1,*i2-1, *val) ;  }
extern "C" void lspecsincrement_(int* i1,int* i2)   { lspecs->increment(*i1-1,*i2-1) ;  }

Mudiar<float> *wgtpairst ;
extern "C" void wgtpairstcreate_(int* n1)      { wgtpairst = new Mudiar<float>(*n1) ;  }
extern "C" void wgtpairstdestroy_(void)                { delete wgtpairst ; }
extern "C" void wgtpairstset_(int* i1,float* val){ wgtpairst->set(*i1-1, *val) ;  }
extern "C" void wgtpairstget_(int* i1,float* val){ wgtpairst->get(*i1-1, *val) ;  }

Mudiar<int> *idpairst ;
extern "C" void idpairstcreate_(int* n1,int* n2)      { idpairst = new Mudiar<int>(*n1,*n2) ;  }
extern "C" void idpairstdestroy_(void)                        { delete idpairst ; }
extern "C" void idpairstset_(int* i1,int* i2,int* val){ idpairst->set(*i1-1,*i2-1, *val) ;  }
extern "C" void idpairstget_(int* i1,int* i2,int* val){ idpairst->get(*i1-1,*i2-1, *val) ;  }

Mudiar<int> *lkfok ;
extern "C" void lkfokcreate_(int* n1, int* n2, int* n3, int* n4, int* n5)       { lkfok = new Mudiar<int>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void lkfokdestroy_(void)                                             { delete lkfok ; }
extern "C" void lkfokset_(int* i1, int* i2, int* i3, int* i4, int* i5, int* val){ lkfok->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void lkfokget_(int* i1, int* i2, int* i3, int* i4, int* i5, int* val){ lkfok->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void lkfokincrement_(int* i1, int* i2, int* i3, int* i4, int* i5)    { lkfok->increment(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1) ; }

//-------------------------------------------------------------------------------------------

Mudiar<float> *emuc ;
extern "C" void createemuc_(int* n1, int* n2, int* n3, int* n4, int* n5)          { emuc = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyemuc_(void)                                                { delete emuc ; }
extern "C" void emucset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuc->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void emucget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuc->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
 
//-------------------------------------------------------------------------------------------

Mudiar<float> *velio ;
extern "C" void createvelio_(int* n1, int* n2, int* n3, int* n4, int* n5)          { velio = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyvelio_(void)                                                { delete velio ; }
extern "C" void velioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { velio->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void velioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { velio->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

Mudiar<float> *bario ;
extern "C" void createbario_(int* n1, int* n2, int* n3, int* n4, int* n5)          { bario = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroybario_(void)                                                { delete bario ; }
extern "C" void barioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { bario->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void barioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { bario->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

Mudiar<float> *epsio ;
extern "C" void createepsio_(int* n1, int* n2, int* n3, int* n4, int* n5)          { epsio = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyepsio_(void)                                                { delete epsio ; }
extern "C" void epsioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { epsio->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void epsioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { epsio->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

Mudiar<float> *emuio ;
extern "C" void createemuio_(int* n1, int* n2, int* n3, int* n4, int* n5)          { emuio = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyemuio_(void)                                                { delete emuio ; }
extern "C" void emuioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuio->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void emuioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuio->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

//-------------------------------------------------------------------------------------------

CCCptl ** cccptl ;

extern "C" void createcccptl_(int* n)
                        {
                          cccptl = new CCCptl* [*n] ;   // here, array of pointers first
                          for(int i=0;i<*n;i++){
                            cccptl[i] = new CCCptl ;
                          }
                        }

extern "C" void destroycccptl_(int* n)
                        {
                          for(int i=0;i<*n;i++){
                            delete cccptl[i] ;
                          }
                          delete [] cccptl ;     // here, array of pointers is deleted last
                        } 

extern "C" void cccptldump_( int* id, int* ist, int* ity, int* ior, int* jor, int* ifr1, int* ifr2
                              , float* tiv1, float* tiv2, float* p1, float* p2, float* p3, float*  p4, float* p5
                              , float* xor1, float* xor2, float* xor3, float* xor4, float* rad
                              , float* des, float* dez, float* qsq, float* zpa1, float* zpa2, float* rin
                              , int* ib1, int* ib2, int* ib3, int* ib4, int* iaa, int* its
                              , int* i )
                        {
                          cccptl[*i]->dump(    *id, *ist, *ity, *ior, *jor, *ifr1, *ifr2
                          , *tiv1, *tiv2, *p1, *p2, *p3, *p4, *p5
                          , *xor1, *xor2, *xor3, *xor4, *rad
                          , *des, *dez, *qsq, *zpa1, *zpa2, *rin
                          , *ib1, *ib2, *ib3, *ib4, *iaa, *its) ;
                           //if(*id==31||*id==25)cout << "ERROR DUMP weird id "<<*i<<"  "<<*id<<"  "<<*ist<<"  "<<*ity<<endl;
                        }
extern "C" void cccptlrestore_( int* id, int* ist, int* ity, int* ior, int* jor, int* ifr1, int* ifr2
                              , float* tiv1, float* tiv2, float* p1, float* p2, float* p3, float*  p4, float* p5
                              , float* xor1, float* xor2, float* xor3, float* xor4, float* rad
                              , float* des, float* dez, float* qsq, float* zpa1, float* zpa2, float* rin
                              , int* ib1, int* ib2, int* ib3, int* ib4, int* iaa, int* its
                              , int* i )
                        {
                          cccptl[*i]->restore( *id, *ist, *ity, *ior, *jor, *ifr1, *ifr2
                          , *tiv1, *tiv2, *p1, *p2, *p3, *p4, *p5
                          , *xor1, *xor2, *xor3, *xor4, *rad
                          , *des, *dez, *qsq, *zpa1, *zpa2, *rin
                          , *ib1, *ib2, *ib3, *ib4, *iaa, *its) ;
                          //if(*id==31||*id==25)cout << "ERROR REST weird id "<<*i<<"  "<<*id<<"  "<<*ist<<"  "<<*ity<<endl;
                        }
extern "C" void  idptlset_(int * i, int * ival)                                                 { cccptl[*i]->setID(*ival); }          
extern "C" void istptlset_(int * i, int * ival)                                                 { cccptl[*i]->setIST(*ival); }          
extern "C" void ityptlset_(int * i, int * ival)                                                 { cccptl[*i]->setITY(*ival); }
extern "C" void iorptlset_(int * i, int * ival)                                                 { cccptl[*i]->setIOR(*ival); }     
extern "C" void jorptlset_(int * i, int * ival)                                                 { cccptl[*i]->setJOR(*ival); }     
extern "C" void ifrptlset_(int * i, int * i1, int * i2)                                         { cccptl[*i]->setIFR(*i1,*i2); }     
extern "C" void desptlset_(int * i, float * val)                                                { cccptl[*i]->setDES(*val); }
extern "C" void radptlset_(int * i, float * val)                                                { cccptl[*i]->setRAD(*val); }
extern "C" void rinptlset_(int * i, float * val)                                                { cccptl[*i]->setRIN(*val); }
extern "C" void   pptlset_(int * i, float * p1, float * p2, float * p3, float * p4, float * p5) { cccptl[*i]->setP(*p1,*p2,*p3,*p4,*p5);}
extern "C" void xorptlset_(int * i, float * x1, float * x2, float * x3, float * x4)             { cccptl[*i]->setXOR(*x1,*x2,*x3,*x4);}
extern "C" void tivptlset_(int * i, float * t1, float * t2)                                     { cccptl[*i]->setTIV(*t1,*t2);}
extern "C" void zpaptlset_(int * i, float * t1, float * t2)                                     { cccptl[*i]->setZPA(*t1,*t2);}
extern "C" void  ibptlset_(int * i, int * i1, int * i2, int * i3, int * i4)                     { cccptl[*i]->setIB(*i1,*i2,*i3,*i4); }     
extern "C" void ibptlset2_(int * i, int * j, int * ival)                                        { cccptl[*i]->set2IB(*j,*ival); }     
 
extern "C" void  idptlget_(int * i, int * ival)                                                 { *ival=cccptl[*i]->getID(); }          
extern "C" void istptlget_(int * i, int * ival)                                                 { *ival=cccptl[*i]->getIST(); }          
extern "C" void ityptlget_(int * i, int * ival)                                                 { *ival=cccptl[*i]->getITY(); }        
extern "C" void iorptlget_(int * i, int * ival)                                                 { *ival=cccptl[*i]->getIOR(); }     
extern "C" void jorptlget_(int * i, int * ival)                                                 { *ival=cccptl[*i]->getJOR(); }     
extern "C" void ifrptlget_(int * i, int * i1, int * i2)                                         { cccptl[*i]->getIFR(*i1,*i2); }     
extern "C" void desptlget_(int * i, float * val)                                                { *val=cccptl[*i]->getDES(); }          
extern "C" void radptlget_(int * i, float * val)                                                { *val=cccptl[*i]->getRAD(); }          
extern "C" void rinptlget_(int * i, float * val)                                                { *val=cccptl[*i]->getRIN(); }          
extern "C" void   pptlget_(int * i, float * p1, float * p2, float * p3, float * p4, float * p5) { cccptl[*i]->getP(*p1,*p2,*p3,*p4,*p5);}
extern "C" void xorptlget_(int * i, float * x1, float * x2, float * x3, float * x4)             { cccptl[*i]->getXOR(*x1,*x2,*x3,*x4);}
extern "C" void tivptlget_(int * i, float * t1, float * t2)                                     { cccptl[*i]->getTIV(*t1,*t2);}
extern "C" void zpaptlget_(int * i, float * t1, float * t2)                                     { cccptl[*i]->getZPA(*t1,*t2);}
extern "C" void  ibptlget_(int * i, int * i1, int * i2, int * i3, int * i4)                     { cccptl[*i]->getIB(*i1,*i2,*i3,*i4); }     
extern "C" void ibptlget2_(int * i, int * j, int * ival)                                        { cccptl[*i]->get2IB(*j,*ival); }     

//-------------------------------------------------------------------------------------------

  /*   OmTab removed after version 3238  */


