#include "ParticlePairsCentralityTrigger.hpp"
#include <cmath>
#include <exception>
#include <string>
#include <vector>

//functions to obtain event and particle information from EPOS at the end of an event simulation
void getNptl(int *numberOfParticles);
void getJorptl(int *indice, int *currentJor);
void getIdptl(int *indice, int *currentId);
void getIstptl(int *indice, int *currentIst);
void getPptl(int *indice, float *px, float *py, float *pz, float *energy, float *mass);
void getBim(float *bim);

/**
 * The aim of this method is the creation of a "histogram" which amounts to counting
 * the number of particle pairs per bin with respect to some variable (observable)
 * defined between some minimum value (minvalue) and some maximum value (maxvalue).
 * The counting is realized using binCounts[binNb][j], where binNb is the bin
 * index, and j an integer between 0 and NR_OF_COLUMNS. binCounts[binNb][0] is
 * automatically filled with the value of the variable "observable" representing
 * the middle of the bin with number binNb. In the present module, one uses the
 * simple choice of binCounts[binNb][j] for j=1 and j=2 being equal, and both
 * counting particles per interval.
 * The function analyze is called at the end of an event simulation. 
 */
void ParticlePairsCentralityTrigger::analyze() {

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
    int selectedId = subanalysis.getSubparametersValue(0);
    float centMin = 0, centMax = 0, massMin = 0, massMax = 0, cent = 0;
    int  selectedIst = 0;
    bool triggerMass = false;
    bool triggerCent = false;
    bool trigger = true;
    if (meaning == "TrgMassTrgCent") { //an example of how to use meaning
      triggerMass = true;
      triggerCent = true;
      selectedIst = subanalysis.getSubparametersValue(1);
      massMin = subanalysis.getSubparametersValue(2);
      massMax = subanalysis.getSubparametersValue(3);
      centMin = subanalysis.getSubparametersValue(4);
      centMax = subanalysis.getSubparametersValue(5);
    } //here one may add other trigger choices

    //get centrality variable, here impact parameter, from EPOS
    getBim(&cent);

    //centrality trigger
    trigger = true;
    if (triggerCent) trigger = (cent >= centMin) and (cent <= centMax);

    if (trigger) {

      //declarations
      int numberOfParticles, currentId, currentIdx, currentIst;
      float px, py, pz, energy, mass, zero=0., value, ptr, mtr, rap; //eta,mom
      float px1, py1, pz1, energy1, mass1;
      float px2, py2, pz2, energy2, mass2;
      //get number of particles
      getNptl(&numberOfParticles);
      //loop over the output particles
      for (int i = 1; i <= numberOfParticles; i++) {
        //get the particle id and particle status
        getIdptl(&i, &currentId);
        getIstptl(&i, &currentIst);
        //select the particles 
        if ((currentId == selectedId || currentId == -1*selectedId) && currentIst == selectedIst) {
          currentIdx=currentId;
          //loop over the output particles
          for (int j = 1; j < i; j++) {
            getIdptl(&j, &currentId);       
            getIstptl(&j, &currentIst);
            if (currentId == -1*currentIdx && currentIst == selectedIst) {
              //get the particles momentum vector, energy, mass
              getPptl(&i, &px1, &py1, &pz1, &energy1, &mass1);
              getPptl(&j, &px2, &py2, &pz2, &energy2, &mass2);
              px=px1+px2;
              py=py1+py2;
              pz=pz1+pz2;
              energy=energy1+energy2;
              //compute transverse momentum, mass, momentum, transverse mass, rapidity, pseudorapidity eta
              ptr = sqrt(pow(px, 2) + pow(py, 2));
              mass=(energy-pz)*(energy+pz)-ptr*ptr;
              mass=sqrt(std::max(zero,mass));
              //mom = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
              mtr = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
              rap = copysign(1., pz) * log((energy + fabs(pz)) / mtr);
              //eta = copysign(1., pz) * log((mom + fabs(pz)) / ptr);
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
                if (triggerMass)
                  trigger = (mass >= massMin) && (mass <= massMax);
                if (trigger) {
                  //get the bin index corresponding to the value value
                  unsigned int binNb =
                      ((value - this->getMinvalue()) / (this->getMaxvalue() - this->getMinvalue())) * this->getNumberOfBins();
                  //increment the corresponding histogram bin
                  binCounts[binNb][1] += 1.;
                  binCounts[binNb][2] += 1.;
                }  
              }
            }
          }//J
        }
      } //i
      
      //consider updated bincounts
      subanalysis.incrementNumberOfTriggeredEvents();//count valid events
      subanalysis.setBinCounts(binCounts);           //put updated binCounts back into subanalysis

    } //centrality trigger

  } //subanalysis
};
