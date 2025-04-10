#include "FlowScalarProductCentralityTrigger.hpp"
#include <cmath>
#include <exception>
#include <string>
#include <vector>
//#include <iostream>
//#include <iomanip>
//using std::setw;

//functions to obtain event and particle information from EPOS at the end of an event simulation
void getNptl(int *numberOfParticles);
void getIdptl(int *index, int *currentId);
void getIstptl(int *index, int *currentIst);
void getPptl(int *index, float *px, float *py, float *pz, float *energy, float *mass);
void getBim(float *bim);
void getCharge(int *index, float *charge);
void getIhadron(int *index, int *ihadron);
void getPolar(float *px, float *py, float *phi);

/**
 * The aim of this method is the creation of a "histogram" which amounts to counting
 * elements necessary to compute flow harmonics based on the scalar product method,
 * per bin with respect to some variable (observable), defined between some minimum 
 * value (minvalue) and some maximum value (maxvalue).
 * The counting is realized using binCounts[binNb][j], where binNb is the bin
 * index, and j an integer between 0 and NR_OF_COLUMNS. binCounts[binNb][0] is
 * automatically filled with the value of the variable "observable" representing
 * the middle of the bin with number binNb. 
 * The function analyze is called at the end of an event simulation. 
 */
void FlowScalarProductCentralityTrigger::analyze() {

  //get information from optns
  std::string observable = this->getObservable();
  float minvalue = this->getMinvalue();
  float maxvalue = this->getMaxvalue();
  int   nbinstot = this->getNumberOfBins();
  int   nbins    = nbinstot - 4 ;
  
  //loop over Subanalyses (defined via subparameters in optns)
  for (Subanalysis &subanalysis : this->getSubanalyses()) {

    std::vector<float> countx(nbins,0.0) ;
    std::vector<float> county(nbins,0.0) ;

    //get the current binCounts from subanalysis
    std::vector<std::array<float, NR_OF_COLUMNS>> binCounts = subanalysis.getBinCounts();

    //get meaning
    std::string meaning = subanalysis.getMeaning();

    //get and interprete subparameters from optns
    float centMin = 0, centMax = 0, cent = 0, phi = 0, qAqB, qAqC, qBqC;
    float etaMinA = 0, etaMaxA = 0, etaMinB = 0, etaMaxB = 0, rapMinB = 0, rapMaxB = 0, etaMinC = 0, etaMaxC = 0;
    int  istMin = 0, istMax = 1, selectedId = 0, iflow = 0;
    bool triggerCent = false;
    bool triggerB = false;    
    bool trigger = true;
    if (meaning == "Id-FlowIndex-A-B1-B2-C-TrgCent") {
      triggerCent = true;
      selectedId = subanalysis.getSubparametersValue(0);
      iflow   = subanalysis.getSubparametersValue(1); //Flow index
      etaMinA = subanalysis.getSubparametersValue(2);  //Range A for q_A
      etaMaxA = subanalysis.getSubparametersValue(3);
      etaMinB = subanalysis.getSubparametersValue(4);  //Range B1 for q_B
      etaMaxB = subanalysis.getSubparametersValue(5);
      rapMinB = subanalysis.getSubparametersValue(6);  //Range B2 for u_B
      rapMaxB = subanalysis.getSubparametersValue(7);
      etaMinC = subanalysis.getSubparametersValue(8);  //Range C for q_C
      etaMaxC = subanalysis.getSubparametersValue(9);
      centMin = subanalysis.getSubparametersValue(10);
      centMax = subanalysis.getSubparametersValue(11);
    } 
    else
    { std::cout << "ERROR FlowScalarProductCentralityTrigger unknown meaning" << std::endl;  exit(1) ;}

    //get centrality variable, here impact parameter, from EPOS
    getBim(&cent);

    //centrality trigger
    trigger = true;
    if (triggerCent) trigger = (cent >= centMin) and (cent <= centMax);

    if (trigger) {

      //declarations
      int numberOfParticles, currentId, currentIst, ihadron;
      float px, py, pz, energy, mass, obs, ptr, mtr, rap, charge, eta, mom;

      //get number of particles
      getNptl(&numberOfParticles);
 
      //First loop over the output particles, compute q-vectors 
      float qvecA1=0, qvecA2=0, countA=0, qvecB1=0, qvecB2=0, countB=0, qvecC1=0, qvecC2=0, countC=0;
      for (int i = 1; i <= numberOfParticles; i++) {
        //get the particle id and particle status
        getIdptl(&i, &currentId);
        getIstptl(&i, &currentIst);
        //select the particles with correct status
        if (currentIst >= istMin && currentIst <= istMax) {
          //get the particle momentum vector, energy, mass, charge, ihadron (=1 if hadron)
          getPptl(&i, &px, &py, &pz, &energy, &mass);
          getCharge(&currentId, &charge);
          getIhadron(&currentId, &ihadron);
          //compute transverse momentum, momentum, transverse mass, rapidity, pseudorapidity eta
          ptr = sqrt(pow(px, 2) + pow(py, 2));
          mom = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
          mtr = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
          rap = copysign(1., pz) * log((energy + fabs(pz)) / mtr);
          eta = copysign(1., pz) * log((mom + fabs(pz)) / ptr);
          getPolar( &px, &py, &phi);
          //Range A
          if( ( eta >= etaMinA ) && ( eta <= etaMaxA ) && (fabs(charge) > 0.1) && (ihadron == 1) && (currentIst == 0) ) {
            qvecA1 = qvecA1 + cos(iflow*phi);
            qvecA2 = qvecA2 + sin(iflow*phi);
            countA = countA + 1;
          }
          //Range B1
          if( ( eta >= etaMinB) && ( eta <= etaMaxB ) && (fabs(charge) > 0.1) && (ihadron == 1) && (currentIst == 0) ) { 
            qvecB1 = qvecB1 + cos(iflow*phi);
            qvecB2 = qvecB2 + sin(iflow*phi);
            countB = countB + 1;
            //std::cout <<"TESTSP/XA"<< setw(12)<<i<< setw(12)<<phi<< setw(12)<<cos(iflow*phi)<< setw(12)<<qvecB1<< setw(12)<<countB <<std::endl;
          }
          //Range C
          if( ( eta >= etaMinC) && ( eta <= etaMaxC ) && (fabs(charge) > 0.1) && (ihadron == 1) && (currentIst == 0) ) {
            qvecC1 = qvecC1 + cos(iflow*phi);
            qvecC2 = qvecC2 + sin(iflow*phi);
            countC = countC + 1;
          }
        }
      } //i
      //normalization
      if(countA>0) qvecA1 = qvecA1 / countA;
      if(countA>0) qvecA2 = qvecA2 / countA;
      if(countB>0) qvecB1 = qvecB1 / countB;
      if(countB>0) qvecB2 = qvecB2 / countB;
      if(countC>0) qvecC1 = qvecC1 / countC;
      if(countC>0) qvecC2 = qvecC2 / countC;

      //products qA*qB etc
      qAqB= qvecA1 * qvecB1 + qvecA2 * qvecB2;
      qAqC= qvecA1 * qvecC1 + qvecA2 * qvecC2;
      qBqC= qvecB1 * qvecC1 + qvecB2 * qvecC2;

      //std::cout <<"TESTSP/XA"<< setw(12)<< qvecA1<< setw(12)<< qvecA2 << setw(12)<< qvecB1 << setw(12)<< qvecB2 << setw(20)<< qAqB << std::endl; 
      //std::cout <<"TESTSP/XA"<< qAqB  <<" "<<  qAqC   <<" "<<  qBqC  <<std::endl;
       
      //Second loop over the output particles, compute u-vectors
      for (int i = 1; i <= numberOfParticles; i++) {
        //get the particle id and particle status
        getIdptl(&i, &currentId);
        getIstptl(&i, &currentIst);
        //select the particles with correct status
        if (currentIst >= istMin && currentIst <= istMax) {
          //get the particle momentum vector, energy, mass
          getPptl(&i, &px, &py, &pz, &energy, &mass);
          //compute transverse momentum, momentum, transverse mass, rapidity, pseudorapidity eta
          ptr = sqrt(pow(px, 2) + pow(py, 2));
          mom = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
          mtr = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
          rap = copysign(1., pz) * log((energy + fabs(pz)) / mtr);
          eta = copysign(1., pz) * log((mom + fabs(pz)) / ptr);
          getPolar( &px, &py, &phi);
          if (observable == "ptr")
            obs = ptr;
          else if (observable == "rapidity")
            obs = rap;
          else
            obs = 1e30;
          int binNb = -1;
          if (obs >= minvalue) binNb = ((obs - minvalue) / (maxvalue - minvalue)) * nbinstot ; //wrong formula for log binning ... need getBinNb   
          if ( (binNb >= 0) && (binNb < nbins) ) {
            triggerB = ( eta >= etaMinB) && ( eta <= etaMaxB ) && ( fabs(currentId) == selectedId ); //particles of interest in range B1
            if(triggerB && rapMinB != rapMaxB) triggerB = ( rap >= rapMinB) && ( rap <= rapMaxB );   //   and in B2 
            if (triggerB ) { 
              //std::cout <<"TESTSP/XA"<< setw(12)<< i<< setw(12)<< currentId << setw(12) << cos(iflow*phi) * qvecC1 + sin(iflow*phi) * qvecC2 << std::endl; 
              countx[binNb] +=  cos(iflow*phi) * qvecC1 + sin(iflow*phi) * qvecC2    ; // uB*qC
              county[binNb] += 1.;
            }
          }
        }
      } //i

      //fill first nbins slots
      for (int k=0; k < nbins; k++){ 
        binCounts[k][1] += countx[k];
        binCounts[k][2] += county[k];   
      }
      //fill extension
      binCounts[nbins][1]   +=  qAqB ;
      binCounts[nbins+1][1] +=  qAqC ;
      binCounts[nbins+2][1] +=  qBqC ;
      binCounts[nbins+3][1] +=  1 ;
      
      //consider updated bincounts
      subanalysis.incrementNumberOfTriggeredEvents();//count valid events
      subanalysis.setBinCounts(binCounts);           //put updated binCounts back into subanalysis

    } //centrality trigger

  } //subanalysis
};
