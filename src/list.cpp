//
// This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
// (See COPYING file for the text of the licence)
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using std::setw;
using std::string;
using std::cout;

void listParticlesFout(int& n,std::ostream& fout);
void getOutputListScreen(int *);
void getOutputListFile(int *, string&);

void getNptl   (int *);
void getIorptl (int *, int *);
void getJorptl (int *, int *);
void getIfrptl (int *, int *, int *);
void getIdptl  (int *, int *);
void getIstptl (int *, int *);
void getItyptl (int *, int *);
void getPptl   (int *, float *, float *, float *, float *, float *);
//void getXorptl (int *, float *, float *, float *, float *);
void getIdpdg  (int *, int *);



/**
 * Output via printing the particle list to the screen and / or into some file, at the end of an event simulation.
 *
 * It is activated in the optns file via the following (one or both) commands:  
 *
 *     print to screen
 *
 *     print into ./MyFile.txt
 *
 * where the name of the file (here ./MyFile.txt) can be chosen freely.
 *
 */

void listParticles(int& n){

   int flag;
   string file;

   getOutputListScreen(&flag);
   if(flag!=0){
     listParticlesFout(n,std::cout);
   }

   getOutputListFile(&flag,file);
   if(flag!=0){
     std::ofstream fout;
     if(n==1)fout.open (file);
     else    fout.open (file, std::ios::app);
     listParticlesFout(n,fout);
   }
}

/**
 * Output via printing the particle list to fout (stdout or file), at the end of an event simulation.
 */

void listParticlesFout(int& n, std::ostream& fout){

   int nptl,ior,jor,ix,ifr1,ifr2,id,idPDG,ist,ity;
   float p1,p2,p3,energy,mass;
   //float xor1,xor2,xor3,time;
   getNptl(&nptl);

   fout<<std::endl;
   fout<<" "<<std::string(112,'#') << std::endl;
   fout<<" "<<"Event number: "<< n <<std::endl;
   fout<<" "<<"Number of particles: "<<nptl<<std::endl;
   fout<<std::endl;

   for (ix=1;ix<=nptl;ix++){
      getIstptl(&ix, &ist);
      //if(ist==25){ continue; } //ist=25 hard partons are special (not final state particles)
      getIorptl (&ix, &ior);
      getJorptl (&ix, &jor);
      if(ist==20||ist==21)jor=0;  //ist=20 or 21 partons have jor=0; jor value here has special meaning   
      getIfrptl (&ix, &ifr1, &ifr2);
      getIdptl  (&ix, &id);
      getIdpdg  (&id, &idPDG);
      getItyptl (&ix, &ity);
      getPptl   (&ix, &p1, &p2, &p3, &energy, &mass);
      //getXorptl (&ix, &xor1, &xor2, &xor3, &time);
      fout<<"     "
        << setw(6) << ior
        << setw(6) << jor
        << setw(6) << ix
        << setw(6) << ifr1
        << setw(6) << ifr2
        << setw(10) << id
        << setw(8) << idPDG
        << setw(4) << ist
        << setw(4) << ity
      //<< setw(12) << p1
      //<< setw(12) << p2
      //<< setw(12) << p3
        << setw(12) << energy
      //<< setw(12) << mass
      //<< setw(12) << xor1
      //<< setw(12) << xor2
      //<< setw(12) << xor3
      //<< setw(12) << time
        << std::endl;
   }
}



/**
 * Output via printing the particle list into the check file, several times,
 * at different stages in the program flow, for each event.
 * Activated via
 *
 *     print * 2
 *
 * in the optns file.
 */

void listParticlesCheck(string text, int i, int n, int m, string file){

    int nptl,ior,jor,ix,ifr1,ifr2,id,idPDG,ist,ity;
    float p1,p2,p3,energy,mass;
    //float xor1,xor2,xor3,time;
    getNptl(&nptl);
    std::ofstream check;
    check.open (file,std::ios::app);
    if(i > 0) {
      check<<std::endl;
      check<<" "<<std::string(74,'#') << std::endl;
      check<<" "<<std::string(20,'#') <<" "<< text<<" "<<std::string(52-i,'#')<<std::endl;
      check<<" "<<std::string(74,'#') << std::endl;
      //check<<" "<<"Number of particles: "<<nptl<<std::endl;
      check<<std::endl;
    }
    check.close ();
    
    if(m==0)return;

    check.open (file,std::ios::app);
    for (ix=n;ix<=m;ix++){
        getIstptl(&ix, &ist);
        //if(ist!=25){ continue; } //ist=25 hard partons are special (not final state particles)
        getIorptl (&ix, &ior);
        getJorptl (&ix, &jor);
        if(ist==20||ist==21)jor=0;  //value here has special meaning   
        getIfrptl (&ix, &ifr1, &ifr2);
        if(text=="list before fragmentation"&&ist/10==2)ifr1=0; //value here has special meaning   
        getIdptl  (&ix, &id);
        getIdpdg  (&id, &idPDG);
        getItyptl (&ix, &ity);
        getPptl   (&ix, &p1, &p2, &p3, &energy, &mass);
        //getxorptl (&ix, &xor1, &xor2, &xor3, &time);
        check<<"     "
          << setw(6) << ior
          << setw(6) << jor
          << setw(6) << ix
          << setw(6) << ifr1
          << setw(6) << ifr2
          << setw(10) << id
          << setw(8) << idPDG
          << setw(4) << ist
          << setw(4) << ity
        //<< setw(12) << p1
        //<< setw(12) << p2
        //<< setw(12) << p3
          << setw(12) << energy
        //<< setw(12) << mass
        //<< setw(12) << xor1
        //<< setw(12) << xor2
        //<< setw(12) << xor3
        //<< setw(12) << time
          << std::endl;
    }
    check.close ();
}

