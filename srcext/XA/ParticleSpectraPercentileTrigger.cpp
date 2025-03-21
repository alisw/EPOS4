#include "ParticleSpectraPercentileTrigger.hpp"
#include <cmath>
#include <ctime>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
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
void ParticleSpectraPercentileTrigger::analyze() {

  //get information from optns
  std::string observable = this->getObservable();
  float minvalue = this->getMinvalue();
  float maxvalue = this->getMaxvalue();
  std::vector<std::vector<float>> importarray = this->getImportarray();
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
    int selectedId = subanalysis.getSubparametersValue(0);
    float rapMin = 0, rapMax = 0, percentMin = 0, percentMax = 0, eta1 = 0, eta2 = 0, eta3 = 0, eta4 = 0;
    bool triggerRap = false;
    bool triggerPercentile = false;
    bool trigger = true;
    if (meaning == "TrgRapTrgPercentile") {
      triggerRap = true;
      triggerPercentile = true;
      rapMin = subanalysis.getSubparametersValue(1);
      rapMax = subanalysis.getSubparametersValue(2);
      percentMax = 100 - subanalysis.getSubparametersValue(3); //percentile(10) means the multiplicity such that 90%
      percentMin = 100 - subanalysis.getSubparametersValue(4); //of all events have a multiplicity less than this value
      eta1 = subanalysis.getSubparametersValue(5);
      eta2 = subanalysis.getSubparametersValue(6);
      eta3 = subanalysis.getSubparametersValue(7);
      eta4 = subanalysis.getSubparametersValue(8);
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
      if (currentIst == 0 || currentIst == 1) {
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
        if (triggerPercentile && (fabs(charge) > 0.1) && (ihadron == 1)) {
          if ((eta >= eta1) && (eta <= eta2)) {
            multiplicity += 1;
          }
          if ((eta >= eta3) && (eta <= eta4)) {
            multiplicity += 1;
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

    //integral over imported probability distribution \int P(x)dx, x = multiplicity 
    float deltaMult = ((importarray.at(1)).at(0) - (importarray.at(0)).at(0));
    float integralProbaDistr = 0;
    for (auto line : importarray) { integralProbaDistr += line.at(1) * deltaMult; }

    //check integralProbaDistr <= 1
    if (integralProbaDistr > 1.001) {
      std::cout << "####### ERROR 26072024a - integralProbaDistr > 1 #######" << std::endl; 
      exit(1);
    }

    //check percentMax <= 100*integralProbaDistr
    if ((percentMax != 100) and (percentMax > 100 * integralProbaDistr)) {
      std::cout << "####### ERROR 26072024b - percentMax > 100*integralProbaDistr #######" << std::endl;
      exit(1);
    }

    //check minbias trigger
    trigger = (minbias >= multMinbias);

    //check percentile trigger
    if (trigger and triggerPercentile) {

      //search percentileMin
      int i = 0;
      float sum = (importarray.at(i)).at(1) * deltaMult;;
      while (sum * 100 < percentMin) {
        i++;
        sum += (importarray.at(i)).at(1) * deltaMult;
      }
      float x1 = (importarray.at(i)).at(0)-deltaMult/2;
      float x2 = (importarray.at(i)).at(0)+deltaMult/2;
      float h2 = sum;
      float h1 = 0;
      if(i > 0) h1 = sum - (importarray.at(i)).at(1);
      float percentileMin=x1+(0.01*percentMin-h1)/(h2-h1)*(x2-x1);

      //search percentileMax
      float percentileMax;
      if (percentMax == 100) { 
        percentileMax = 1e30; 
      } else {
        while (sum * 100 < percentMax) {
          i++;
          sum += (importarray.at(i)).at(1) * deltaMult;
        }
        x1 = (importarray.at(i)).at(0)-deltaMult/2;
        x2 = (importarray.at(i)).at(0)+deltaMult/2;
        h2 = sum;
        h1 = 0;
        if(i > 0) h1 = sum - (importarray.at(i)).at(1);
        percentileMax=x1+(0.01*percentMax-h1)/(h2-h1)*(x2-x1);
      }
      //std::cout << x1   << "   "<< x2   << "   "<<  h1  << "   "<<  h2  << "   "<<  percentileMax <<   std::endl;
      trigger = (multiplicity >= percentileMin ) and (multiplicity <= percentileMax );
      //std::cout << percentileMin   << "   "<< percentileMax << "   "<< multiplicity  <<std::endl;
    }

    //consider updated bincounts when trigger conditions are satisfied
    if (trigger) {
      subanalysis.incrementNumberOfTriggeredEvents(); //count valid events
      subanalysis.setBinCounts(binCounts);            //put updated binCounts back into subanalysis
    }
  }
};
