#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define Nbin 20
#define pi 3.14159

// declaring global variables

int Nev = 0;   //counts events
float delbin=2*pi/Nbin; 
float distr[Nbin];
void finish_program_4_();

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   the function doz4  activates (or not) the other functions in this file
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int doz4()
{ 
  return 0 ;   /* put "return 1" to activate the following routines */  
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   the function transfer_event_4_
//   is called at the end of each event simulation
//    it transfers the information about all the particles  in the event    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void transfer_event_4_(int *izmode, float *cnt, int *nptevt
                      ,int *id, int *ist, int *ity, int *ior, int *jor
                      ,float *px, float *py, float *pz, float *en, float *am
                      ,float *x, float *y, float *z, float *t)
{

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! All arguments of transfer_event_4_ are pointers.
  ! *izmode, *cnt, *nptevt refer to simple variables,
  ! all others to arrays with particle properties. Array dimension : *nptevt
  ! *izmode is an integer defing the meaning of the centrality variable *cnt :
  !    1 : impact parameter
  !    2 : number of Pomerons
  ! *nptevt is the number of particles in the event.
  ! Meaning of the arrays elements :
  !  id[i] .................... id of particle i
  !  ist[i] ................... status and particle type of particle i
  !  ity[i] ................... type of particle origin of particle i
  !  ior[i]-1 ................... mother of particle i
  !  jor[i]-1 ................... father of particle i
  !  px[i], py[i], pz[i] ...... momentum of particle i
  !  en[i] .................... energy  of particle i
  !  am[i] .................... mass  of particle i
  !  x[i], y[i], z[i], t[i] ... space-time of particle i
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  if(doz4()!=1)return;
  Nev++;

  int i,j,n;

  // initialization

  if(Nev==1)
  { 
      for(n=0;n<Nbin;n++)distr[n]=0;
  }

  // print event info

  //printf("----------- Event number %d-----------------\n",Nev);
  //printf("izmode = %d  cnt = %f  nptevt = %d \n",*izmode,*cnt,*nptevt);

  // correlation
  
  float delphi;
/*
  //loop over first particle:
  for (i = 0; i<*nptevt;i++)
  {
    float pt1=sqrt(px[i]*px[i]+py[i]*py[i]);
    float p1=sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
    //eta is pseudorapidity, find definition via Google 
    float eta1=pz[i]/fabs(pz[i])*log( (p1+fabs(pz[i])) /pt1 );
    //trigger conditions for particle 1:
    if(ist[i]==0 && pt1>0.001 && fabs(eta1)<2)
    {
      //loop over second particle:
      for (j=0; j<i;j++)
      {
        float pt2=sqrt(px[j]*px[j]+py[j]*py[j]);
        float p2=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
        float eta2=pz[j]/fabs(pz[j])*log( (p2+fabs(pz[j])) /pt2 );
        //trigger conditions for particle 2:
        if(ist[j]==0 && pt2>0.001 && fabs(eta2)<2)
        {
              float deleta=((eta2)-(eta1));
              delphi = py[j]/fabs(py[j])*acos((px[j])/pt2) - py[i]/fabs(py[i])*acos((px[i])/pt1)  ;
              if(delphi>pi)delphi=delphi-2*pi;
              if(delphi<-pi)delphi=delphi+2*pi;
              // printf("  %d  %d     %f   %f   %f   \n",c,b,deleta, eta1, eta2 );
              int nbin=(delphi+pi)/delbin; //bin number
              if(nbin>=0&&nbin<Nbin)distr[nbin]++;
        }
      }
    }
  }
*/
     
      int kk,k1,k2,k3,kkc=0;
      for (j=0; j<*nptevt;j++)
      {
         if(abs(id[j]) < 10000 && abs(id[j]) > 99) 
         {  
            kk= floor( (abs(id[j]) % 10000)/10) ;
            k1=kk % 10; 
            k2=(int)floor(kk/10) % 10;
            k3=(int)floor(kk/100) % 10;
            if(k1==4||k2==4||k3==4)
            { 
               kkc++;
               //printf("%d   %d %d %d\n",id[j],k3,k2,k1);
            }
         }
      }

      //loop over second particle:
      if(kkc>0)
      {
       //printf("   ----------- hit !!! ------------- \n");
       for (j=0; j<*nptevt;j++)
       {
        float pt2=sqrt(px[j]*px[j]+py[j]*py[j]);
        float p2=sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
        float eta2=pz[j]/fabs(pz[j])*log( (p2+fabs(pz[j])) /pt2 );
        //trigger conditions for particle 2:
        int ida=abs(id[j]);
        if( (ida==110||ida==120||ida==1120||ida==1220||ida==130||ida==12||ida==14)
        && ist[j]==0  && fabs(eta2)<1)
        {
              delphi = py[j]/fabs(py[j])*acos((px[j])/pt2)  ;
              if(delphi>pi)delphi=delphi-2*pi;
              if(delphi<-pi)delphi=delphi+2*pi;
              // printf("  %d  %d     %f   %f   %f   \n",c,b,deleta, eta1, eta2 );
              int nbin=(delphi+pi)/delbin; //bin number
              if(nbin>=0&&nbin<Nbin)distr[nbin]++;
        }
       }
      }

  // For testing, later comment following line
  //if(Nev==2000){finish_program_4_();exit(2);}
 
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  the function finish_program_4_
//  is called at the very end of the run
//  It can be used to print final distributions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void finish_program_4_()
 {
  int j;
  if(doz4()!=1)return;


  printf("\n---------- distr after %d events ----------\n",Nev);
  for(j=0;j<Nbin;j++)
  {
    printf("  %f %f \n",(j+0.5)*delbin,distr[j]);
  }

}

