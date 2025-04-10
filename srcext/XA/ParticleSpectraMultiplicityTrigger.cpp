#include "ParticleSpectraMultiplicityTrigger.hpp"
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
The aim of this method is the creation of a "histogram" which amounts to counting
the number of particles per bin with respect to some variable (observable)
defined between some minimum value (minvalue) and some maximum value (maxvalue).
The counting is realized using binCounts[binNb][j], where binNb is the bin
index, and j an integer between 0 and NR_OF_COLUMNS. binCounts[binNb][0] is
automatically filled with the value of the variable "observable" representing
the middle of the bin with number binNb. In the present module, one uses the
simple choice of binCounts[binNb][j] for j=1 and j=2 being equal, and both
counting particles per interval.
The function analyze is called at the end of an event simulation. 
 */
void ParticleSpectraMultiplicityTrigger::analyze() {

  //get information from optns
  std::string observable = this->getObservable();
  float minvalue = this->getMinvalue();
  float maxvalue = this->getMaxvalue();
  std::vector<float> moreparam = this->getMoreparameters();
  float eta1Minbias = moreparam[0];
  float eta2Minbias = moreparam[1];
  float eta3Minbias = moreparam[2];  //0 if not
  float eta4Minbias = moreparam[3];  //0 needed
  float multMinbias = moreparam[4];

  //loop over subanalyses (defined via subparameters in optns)
  for (Subanalysis &subanalysis : this->getSubanalyses()) {


    //get the current binCounts from subanalysis
    std::vector<std::array<float, NR_OF_COLUMNS>> binCounts = subanalysis.getBinCounts();

    //get meaning from optns
    std::string meaning = subanalysis.getMeaning();

    //get and interprete subparameters from optns
    float rapMin = 0, rapMax = 0, ptrMin = 0, ptrMax = 0, multMin = 0, multMax = 0, etaMinMult = 0, etaMaxMult = 0;
    int selectedId = 0, istMin = 0, istMax = 1; 
    bool triggerRap = false;
    bool triggerPtr = false;
    bool triggerMult = false;
    bool trigger = true;
    if (meaning == "TrgRap-TrgMult") { //an example of how to use meaning
      triggerRap = true;
      triggerMult = true;
      selectedId = subanalysis.getSubparametersValue(0);
      rapMin = subanalysis.getSubparametersValue(1);
      rapMax = subanalysis.getSubparametersValue(2);
      multMin = subanalysis.getSubparametersValue(3);
      multMax = subanalysis.getSubparametersValue(4);
      etaMinMult = subanalysis.getSubparametersValue(5);
      etaMaxMult = subanalysis.getSubparametersValue(6);
    } else if (meaning == "TrgMult") {
      triggerMult = true;
      selectedId = subanalysis.getSubparametersValue(0);
      multMin = subanalysis.getSubparametersValue(1);
      multMax = subanalysis.getSubparametersValue(2);
      etaMinMult = subanalysis.getSubparametersValue(3);
      etaMaxMult = subanalysis.getSubparametersValue(4);
    } else if (meaning == "TrgRap") {
      triggerRap = true;
      selectedId = subanalysis.getSubparametersValue(0);
      rapMin = subanalysis.getSubparametersValue(1);
      rapMax = subanalysis.getSubparametersValue(2);
    } else if (meaning == "TrgPtr") {
      triggerPtr = true;
      selectedId = subanalysis.getSubparametersValue(0);
      ptrMin = subanalysis.getSubparametersValue(1);
      ptrMax = subanalysis.getSubparametersValue(2);
    } else if (meaning == "TrgIst-TrgRap") {
      triggerRap = true;
      selectedId = subanalysis.getSubparametersValue(0);
      istMin  = subanalysis.getSubparametersValue(1);
      istMax  = subanalysis.getSubparametersValue(2);
      rapMin = subanalysis.getSubparametersValue(3);
      rapMax = subanalysis.getSubparametersValue(4);
    } else if (meaning == "TrgIst-TrgPtr") {
      triggerPtr = true;
      selectedId = subanalysis.getSubparametersValue(0);
      istMin  = subanalysis.getSubparametersValue(1);
      istMax  = subanalysis.getSubparametersValue(2);
      ptrMin = subanalysis.getSubparametersValue(3);
      ptrMax = subanalysis.getSubparametersValue(4);
    }

    //declarations
    int numberOfParticles, currentId, currentIst, ihadron;
    float px, py, pz, energy, mass, mom, charge, obs, ptr, mtr, rap, eta;
    unsigned int multiplicity = 0;
    unsigned int minbias = 0;

    //get number of particles
    getNptl(&numberOfParticles);
    //loop over the output particles
    for (int i = 1; i <= numberOfParticles; i++) {
      //get the particle id and particle status
      getIdptl(&i, &currentId);
      getIstptl(&i, &currentIst);
      //select the particles with status 0 (last generation)
      if (currentIst >= istMin && currentIst <= istMax) {
        //get the particle momentum vector, energy, mass, charge, ihadron (=1 if hadron)
        getPptl(&i, &px, &py, &pz, &energy, &mass);
        getcharge_(&currentId, &charge);
        getihadron_(&currentId, &ihadron);
        //compute transverse momentum, momentum, transverse mass, rapidity, pseudorapidity eta
        ptr = sqrt(pow(px, 2) + pow(py, 2));
        mom = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
        mtr = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
        rap = copysign(1., pz) * log((energy + fabs(pz)) / mtr);
        eta = copysign(1., pz) * log((mom + fabs(pz)) / ptr);
        if (currentId == selectedId) {
          if (observable == "ptr")
            obs = ptr;
          else if (observable == "rapidity")
            obs = rap;
          else
            obs = 1e30;
          //select the particles with obs in range [minvalue,maxvalue]
          if ((obs >= minvalue) && (obs <= maxvalue)) {
            //check trigger conditions
            trigger = true;
            if (triggerRap)
              trigger = (rap >= rapMin) && (rap <= rapMax);
            if (triggerPtr)
              trigger = (ptr >= ptrMin) && (ptr <= ptrMax);
            if (trigger) {
              //get the bin index corresponding to the obs value
              unsigned int binNb =
                  ((obs - this->getMinvalue()) / (this->getMaxvalue() - this->getMinvalue())) * this->getNumberOfBins();
              //increment the corresponding histogram bin
              binCounts[binNb][1] += 1.;
              binCounts[binNb][2] += 1.;
            }
          }
        }
        //count particles for charged hadron multiplicity
        if (triggerMult && (fabs(charge) > 0.1) && (ihadron == 1)) {
          if ((eta >= etaMinMult) && (eta <= etaMaxMult)) {
            multiplicity += 1;
          }
        }
        //count particles for minimum bias condition
        if ((fabs(charge) > 0.1) && (ihadron == 1)) {
          if ((eta >= eta1Minbias) && (eta <= eta2Minbias)) { minbias += 1; }
          if(eta3Minbias != eta4Minbias) {
            if ((eta >= eta3Minbias) && (eta <= eta4Minbias)) { minbias += 1; }
          }
        }
      }
    }

    //check min bias trigger
    trigger = (minbias >= multMinbias);

    //check mult trigger
    if (trigger and triggerMult)
      trigger = (multiplicity >= multMin) and (multiplicity <= multMax);

    //count when trigger conditions are satisfied
    if (trigger) {
      subanalysis.incrementNumberOfTriggeredEvents();//count valid events
      subanalysis.setBinCounts(binCounts);           //put updated binCounts back into subanalysis
    }
  }
};
