#include "MultiplicityDistributions.hpp"
#include <cmath>
#include <exception>
#include <string>
#include <vector>

//functions to obtain event and particle information from EPOS at the end of an event simulation
void getNptl(int *numberOfParticles);
void getIdptl(int *indice, int *currentId);
void getIstptl(int *indice, int *currentIst);
void getPptl(int *indice, float *px, float *py, float *pz, float *energy, float *mass);
extern "C" {
void getcharge_(int *currentId, float *charge);
void getihadron_(int *currentId, int *ihadron);
}

/**
 * The aim of this method is the creation of a "histogram" which amounts to counting 
 * the number of events per bin with respect to some variable (observable) representing 
 * the event multiplicity, defined between some minimum value (minvalue) and some 
 * maximum value (maxvalue). The counting is realized using binCounts[binNb][j], where
 * binNb is the bin index, and j an integer between 0 and NR_OF_COLUMNS. binCounts[binNb][0]
 * is automatically filled with the value of the variable "observable" representing the 
 * middle of the bin with number binNb. In the present module, one uses the simple choice
 * of binCounts[binNb][j] for j=1 and j=2 being equal, and both representing counts per interval.
 * The function analyze is called at the end of an event simulation. 
 */
void MultiplicityDistributions::analyze() {

  //get information from optns
  std::string observable = this->getObservable();
  float minvalue = this->getMinvalue();
  float maxvalue = this->getMaxvalue();
  std::vector<float> moreparam = this->getMoreparameters();
  float eta1Minbias = moreparam[0];
  float eta2Minbias = moreparam[1];
  float multMinbias = moreparam[2];

  //loop over subanalyses (defined via subparameters in optns)
  for (Subanalysis &subanalysis : this->getSubanalyses()) {

    //get the current binCounts from subanalysis
    std::vector<std::array<float, NR_OF_COLUMNS>> binCounts = subanalysis.getBinCounts();
    
    //get meaning from optns
    std::string meaning = subanalysis.getMeaning();

    //get and interprete subparameters from optns
    float eta1=0,eta2=0,eta3=0,eta4=0;
    if(meaning == "AliceV0") { 
      eta1 = subanalysis.getSubparametersValue(0);
      eta2 = subanalysis.getSubparametersValue(1);
      eta3 = subanalysis.getSubparametersValue(2);
      eta4 = subanalysis.getSubparametersValue(3);
    }

    int numberOfParticles, currentId, currentIst, ihadron;
    float px, py, pz, energy, mass, mom, charge, obs=0., ptr, eta; //mtr,rap
    unsigned int multiplicity = 0;
    unsigned int minbias = 0;

    //get number of particles
    getNptl(&numberOfParticles);
    //loop over the particles
    for (int i = 1; i <= numberOfParticles; i++) {
      //get the particle id and particle status
      getIdptl(&i, &currentId);
      getIstptl(&i, &currentIst);
      //select the particles with status 0 (last generation)
      if (currentIst == 0) {
        //get the particle momentum vector, energy, mass, charge, ihadron (=1 if hadron)
        getPptl(&i, &px, &py, &pz, &energy, &mass);
        getcharge_(&currentId,&charge);   //charge
        getihadron_(&currentId,&ihadron); //1 if hadron
        //compute transverse momentum, momentum, transverse mass, rapidity, pseudorapidity eta
        ptr = sqrt(pow(px, 2) + pow(py, 2));
        mom = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
      //mtr = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
      //rap = copysign(1., pz) * log((energy + fabs(pz)) / mtr);
        eta = copysign(1., pz) * log((mom + fabs(pz)) / ptr);
        //count particles for charged hadron multiplicity
        if( meaning == "AliceV0" ) {
          if( (fabs(charge)>0.1) && (ihadron==1) ) {
            if ((eta >= eta1) && (eta <= eta2)) { multiplicity += 1; }
            if ((eta >= eta3) && (eta <= eta4)) { multiplicity += 1; }         
          }
        }
        //count particles for minimum bias condition
        if ((fabs(charge) > 0.1) && (ihadron == 1)) {
          if ((eta >= eta1Minbias) && (eta <= eta2Minbias)) {
            minbias += 1;
          }
        }
      }
    }

    //check min bias trigger
    bool trigger = (minbias >= multMinbias);

    //count 
    if( observable == "multiplicity") obs = multiplicity ;
    if ( trigger and (obs >= minvalue) && (obs <= maxvalue)) {
      //get the bin index corresponding to the obs value
      unsigned int binNb = ((obs - this->getMinvalue()) / (this->getMaxvalue() - this->getMinvalue())) * this->getNumberOfBins();      
      //increment the corresponding histogram bin
      binCounts[binNb][1] += 1.;
      binCounts[binNb][2] += 1.;
    }
    
    //consider updated bincounts when trigger conditions are satisfied
    if (trigger) {
      subanalysis.incrementNumberOfTriggeredEvents();  //count valid events
      subanalysis.setBinCounts(binCounts);             //put updated binCounts back into subanalysis
    }
  }
};
