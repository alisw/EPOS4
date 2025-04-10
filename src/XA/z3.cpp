#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>

#include "z3.h"

int imop,ifap,ip,ic1p,ic2p,idp,istp,ityp ;
float pxp, pyp, pzp, ep, xp, yp, zp, tp ;

using namespace std ;

//--------------------------------------------------------------------------------------------
   void putptl_(int *_imop, int *_ifap, int *_ip, int *_ic1p, int *_ic2p
              , int *_idp, int *_istp, int *_ityp
              , float *_pxp, float *_pyp, float *_pzp, float *_ep
              , float *_xp, float *_yp, float *_zp, float *_tp) 
//--------------------------------------------------------------------------------------------
{

imop = *_imop;
ifap = *_ifap;
ip   = *_ip ; 
ic1p = *_ic1p;
ic2p = *_ic2p;
idp  = *_idp ;
istp = *_istp ;
ityp = *_ityp ;
ep   = *_ep ;
if(istp == 21)ifap=0; //this variable contains different information in case of 21

cout <<imop<<" "<<ifap<<"     "<<ip<<"     "<<ic1p<<" "<<ic2p<<"           "
<<idp<<" "<<istp<<" "<<ityp<<"                "<<ep<< endl ;

}

//--------------------------------------------------------------------------------------------
   void endevt_()
//--------------------------------------------------------------------------------------------
{

cout << " =====================> end event "  << endl ;

}

//--------------------------------------------------------------------------------------------
   void endall_()
//--------------------------------------------------------------------------------------------
{

cout << " ================================> end "  << endl ;
ofstream fout ("zzzTest.txt") ;

fout << " final results" << endl; 


}
