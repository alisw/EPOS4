#include "PtSpectraMultTrigger.hpp"
#include <cmath>
#include <exception>
#include <string>
#include <vector>

//functions to obtain event and particle information from EPOS at the end of an event simulation
void getNptl(int *numberOfParticles);
void getIdptl(int *indice, int *currentId);
void getIstptl(int *indice, int *currentIst);
void getPptl(int *indice, float *px, float *py, float *pz, float *energy, float *mass);

/**
 * The aim of this method is the creation of a "histogram" which amounts to counting
 * the number of particles per bin with respect to the transverse momentum
 * defined between some minimum value (ptMin) and some maximum value (ptMax).
 * The counting is realized using binCounts[binNb][j], where binNb is the bin
 * index, and j an integer between 0 and NR_OF_COLUMNS. binCounts[binNb][0] is
 * automatically filled with the value of the variable "observable" representing
 * the middle of the bin with number binNb. In the present module, one uses the
 * simple choice of binCounts[binNb][j] for j=1 and j=2 being equal, and both
 * counting particles per interval.
 * The function analyze is called at the end of an event simulation. 
 */
void PtSpectraMultTrigger::analyze() {

  //get information from optns
  float ptMin = this->getMinvalue();
  float ptMax = this->getMaxvalue();

  //loop over subanalyses (defined via subparameters in optns)
  for (Subanalysis &subanalysis : this->getSubanalyses()) {

    //get the current binCounts from subanalysis
    std::vector<std::array<float, NR_OF_COLUMNS>> binCounts = subanalysis.getBinCounts();

    //get and interprete subparameters from optns
    int selectedId = subanalysis.getSubparametersValue(0);
    float rapidityMin = subanalysis.getSubparametersValue(1);
    float rapidityMax = subanalysis.getSubparametersValue(2);
    float multMin = subanalysis.getSubparametersValue(3);
    float multMax = subanalysis.getSubparametersValue(4);
    float rapidityMinMult = subanalysis.getSubparametersValue(5);
    float rapidityMaxMult = subanalysis.getSubparametersValue(6);

    //declarations
    int numberOfParticles = 0;
    int currentId = 0;
    int currentIst = 0;
    float px, py, pz, energy, mass, transverseMomentum, transverseMass, rapidity;
    unsigned int multiplicity = 0;

    //get number of particles
    getNptl(&numberOfParticles);
    //loop over the output particles
    for (int i = 1; i <= numberOfParticles; i++) {
      //get the particle id and particle status
      getIdptl(&i, &currentId);
      getIstptl(&i, &currentIst);
      //select the particles with status 0 (last generation)
      if (currentIst == 0 || currentIst == 1) {
        //get the particle momentum, energy and mass
        getPptl(&i, &px, &py, &pz, &energy, &mass);
        //compute the transverse momentum, transverse mass, rapidity for this particle
        transverseMomentum = sqrt(pow(px, 2) + pow(py, 2));
        transverseMass = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
        rapidity = copysign(1., pz) * log((energy + fabs(pz)) / transverseMass);
        if (currentId == selectedId) {
          //select the particles with transverse momentum in range [ptMin,ptMax]
          if ((transverseMomentum >= ptMin) && (transverseMomentum <= ptMax)) {
            //select the particles with rapidity in range [rapidityMin,rapidityMax]
            if ((rapidity >= rapidityMin) && (rapidity <= rapidityMax)) {
              //get the bin index corresponding to the transverse momentum value
              unsigned int binNb = this->getBinNb(transverseMomentum);
              //increment the corresponding histogram bin
              binCounts[binNb][1] += 1.;
              binCounts[binNb][2] += 1.;
            }
          }
        }
        //count particles for charged hadron multiplicity
        if ((rapidity >= rapidityMinMult) && (rapidity <= rapidityMaxMult)) {
          multiplicity += 1;
        }
      }
    }
    
    //check multiplicity trigger
    bool trigger = (multiplicity >= multMin) and (multiplicity <= multMax);
    
    //count when trigger conditions are satisfied
    if (trigger) {
      subanalysis.incrementNumberOfTriggeredEvents(); //count valid events
      subanalysis.setBinCounts(binCounts);            //put updated binCounts back into subanalysis
    }
  }
};
