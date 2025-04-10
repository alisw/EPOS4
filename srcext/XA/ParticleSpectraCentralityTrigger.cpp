#include "ParticleSpectraCentralityTrigger.hpp"
#include <cmath>
#include <exception>
#include <string>
#include <vector>

//functions to obtain event and particle information from EPOS at the end of an event simulation
void getNptl(int *numberOfParticles);
void getIdptl(int *indice, int *currentId);
void getIstptl(int *indice, int *currentIst);
void getPptl(int *indice, float *px, float *py, float *pz, float *energy, float *mass);
void getBim(float *bim);

/**
 * The aim of this method is the creation of a "histogram" which amounts to counting
 * the number of particles per bin with respect to some variable (observable)
 * defined between some minimum value (minvalue) and some maximum value (maxvalue).
 * The counting is realized using binCounts[binNb][j], where binNb is the bin
 * index, and j an integer between 0 and NR_OF_COLUMNS. binCounts[binNb][0] is
 * automatically filled with the value of the variable "observable" representing
 * the middle of the bin with number binNb. In the present module, one uses the
 * simple choice of binCounts[binNb][j] for j=1 and j=2 being equal, and both
 * counting particles per interval.
 * The function analyze is called at the end of an event simulation. 
 */
void ParticleSpectraCentralityTrigger::analyze() {

  //get information from optns
  std::string observable = this->getObservable();
  float minvalue = this->getMinvalue();
  float maxvalue = this->getMaxvalue();

  //loop over Subanalyses (defined via subparameters in optns)
  for (Subanalysis &subanalysis : this->getSubanalyses()) {

    //get the current binCounts from subanalysis
    std::vector<std::array<float, NR_OF_COLUMNS>> binCounts = subanalysis.getBinCounts();

    //get meaning
    std::string meaning = subanalysis.getMeaning();

    //get and interprete subparameters from optns
    float centMin = 0, centMax = 0, rapMin = 0, rapMax = 0, cent = 0;
    int  nid = 1, istMin = 0, istMax = 1, numberOfIds = 1, selectedIds[9] = {};
    bool triggerRap = false;
    bool triggerCent = false;
    bool triggerId = false;
    bool trigger = true;
    if (meaning == "TrgRap-TrgCent") {
      triggerRap = true;
      triggerCent = true;
      selectedIds[0] = subanalysis.getSubparametersValue(0);
      rapMin = subanalysis.getSubparametersValue(1);
      rapMax = subanalysis.getSubparametersValue(2);
      centMin = subanalysis.getSubparametersValue(3);
      centMax = subanalysis.getSubparametersValue(4);
    } else if (meaning == "TwoIds-TrgRap-TrgCent") {
      triggerRap = true;
      triggerCent = true;
      selectedIds[0] = subanalysis.getSubparametersValue(0);
      nid = subanalysis.getSubparametersValue(1);
      rapMin = subanalysis.getSubparametersValue(2);
      rapMax = subanalysis.getSubparametersValue(3);
      centMin = subanalysis.getSubparametersValue(4);
      centMax = subanalysis.getSubparametersValue(5);
    } else if (meaning == "MultIds-TrgIst-TrgRap-TrgCent") { // allows specifying several ids 
      triggerRap = true;
      triggerCent = true;
      numberOfIds = subanalysis.getSubparametersValue(0);
      for (int k=0;k<numberOfIds;k++){
        selectedIds[k] = subanalysis.getSubparametersValue(k+1);      
      }
      istMin  = subanalysis.getSubparametersValue(numberOfIds+1);
      istMax  = subanalysis.getSubparametersValue(numberOfIds+2);
      rapMin  = subanalysis.getSubparametersValue(numberOfIds+3);
      rapMax  = subanalysis.getSubparametersValue(numberOfIds+4);
      centMin = subanalysis.getSubparametersValue(numberOfIds+5);
      centMax = subanalysis.getSubparametersValue(numberOfIds+6);
    } //here one may add other trigger choices

    //get centrality variable, here impact parameter, from EPOS
    getBim(&cent);

    //centrality trigger
    trigger = true;
    if (triggerCent) trigger = (cent >= centMin) and (cent <= centMax);

    if (trigger) {

      //declarations
      int numberOfParticles, currentId, currentIst;
      float px, py, pz, energy, mass, value, ptr, mtr, rap; //eta,mom

      //get number of particles
      getNptl(&numberOfParticles);
      //loop over the output particles
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
          //mom = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
          mtr = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
          rap = copysign(1., pz) * log((energy + fabs(pz)) / mtr);
          //eta = copysign(1., pz) * log((mom + fabs(pz)) / ptr);
          triggerId = false;
          for (int k=0;k<numberOfIds;k++){
            if (currentId == selectedIds[k]) triggerId = true;
            if (nid == 2 && currentId == -1*selectedIds[k]) triggerId = true;
          }
          if (triggerId) {
            if (observable == "ptr")
              value = ptr;
            else if (observable == "rapidity")
              value = rap;
            else
              value = 1e30;
            //select the particles with value in range [minvalue,maxvalue]
            if ((value >= minvalue) && (value <= maxvalue)) {
              //check trigger conditions
              trigger = true;
              if (triggerRap)
                trigger = (rap >= rapMin) && (rap <= rapMax);
              if (trigger) {
                //get the bin index corresponding to the value value
                unsigned int binNb =
                    ((value - this->getMinvalue()) / (this->getMaxvalue() - this->getMinvalue())) * this->getNumberOfBins();
                //increment the corresponding histogram bin
                binCounts[binNb][1] += 1./nid;
                binCounts[binNb][2] += 1./nid;
              }
            }
          }
        }
      } //i
      
      //consider updated bincounts
      subanalysis.incrementNumberOfTriggeredEvents();//count valid events
      subanalysis.setBinCounts(binCounts);           //put updated binCounts back into subanalysis

    } //centrality trigger

  } //subanalysis
};
