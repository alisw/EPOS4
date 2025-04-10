#include "MultiplicityDependences.hpp"
#include <cmath>
#include <exception>
#include <string>
#include <vector>
//#include <iomanip>
//using std::setw;

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
void MultiplicityDependences::analyze() {

  //get information from optns
  std::string observable = this->getObservable();
  float minvalue = this->getMinvalue();
  float maxvalue = this->getMaxvalue();
  int   nbinstot = this->getNumberOfBins();
  int   nbins    = nbinstot - 3 ;
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
    int selectedId = 0, istMin = 0, istMax = 1;
    float rap1Particle = 0, rap2Particle = 0, eta1Multiplicity = 0, eta2Multiplicity = 0;
    if(meaning == "Id-Ist12-Rap12Particle-Eta12Multiplicity") {
      selectedId = subanalysis.getSubparametersValue(0);
      istMin  = subanalysis.getSubparametersValue(1);
      istMax  = subanalysis.getSubparametersValue(2);
      rap1Particle = subanalysis.getSubparametersValue(3);
      rap2Particle = subanalysis.getSubparametersValue(4);
      eta1Multiplicity = subanalysis.getSubparametersValue(5);
      eta2Multiplicity = subanalysis.getSubparametersValue(6);
    }
    else
    { std::cout << "ERROR MultiplicityDependences unknown meaning" << std::endl;  exit(1) ;}

    int numberOfParticles, currentId, currentIst, ihadron;
    float px, py, pz, energy, mass, mom, charge, obs=0., ptr, eta, mtr,rap;
    unsigned int multiplicity = 0;
    unsigned int minbias = 0;
    unsigned int npoi = 0;

    //get number of particles
    getNptl(&numberOfParticles);
    //loop over the particles
    for (int i = 1; i <= numberOfParticles; i++) {
      //get the particle id and particle status
      getIdptl(&i, &currentId);
      getIstptl(&i, &currentIst);
      //get the particle momentum vector, energy, mass, charge, ihadron (=1 if hadron)
      getPptl(&i, &px, &py, &pz, &energy, &mass);
      getcharge_(&currentId,&charge);   //charge
      getihadron_(&currentId,&ihadron); //1 if hadron
      //compute transverse momentum, momentum, transverse mass, rapidity, pseudorapidity eta
      ptr = sqrt(pow(px, 2) + pow(py, 2));
      mom = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
      mtr = sqrt(pow(mass, 2) + pow(px, 2) + pow(py, 2));
      rap = copysign(1., pz) * log((energy + fabs(pz)) / mtr);
      eta = copysign(1., pz) * log((mom + fabs(pz)) / ptr);

      //count particles for charged hadron multiplicity
      if( meaning == "Id-Ist12-Rap12Particle-Eta12Multiplicity" ) {
        if (currentIst == 0 ) {
          if( (fabs(charge)>0.1) && (ihadron==1) ) {
            if ((eta >= eta1Multiplicity) && (eta <= eta2Multiplicity)) { multiplicity += 1; }         
          }
        }
      }
      else
      { std::cout << "ERROR MultiplicityDependences unknown meaning (2)" << std::endl;  exit(1) ;}

      //count particles for minimum bias condition
      if (currentIst == 0 ) {
        if(eta1Minbias != eta2Minbias) {
         if ((fabs(charge) > 0.1) && (ihadron == 1)) {
          if ((eta >= eta1Minbias) && (eta <= eta2Minbias)) { minbias += 1; }
            if(eta3Minbias != eta4Minbias) {
              if ((eta >= eta3Minbias) && (eta <= eta4Minbias)) { minbias += 1; }
            }
          } 
        }
      }

      //count particles of interest (poi)
      if ( currentId == selectedId ) {
        if (currentIst >= istMin && currentIst <= istMax) {
          if ((rap >= rap1Particle) && (rap <= rap2Particle)) {
            npoi += 1;
          }
        } 
      }
    }

    //check min bias trigger
    bool trigger = (minbias >= multMinbias);

    //count 
    if( observable == "multiplicity") obs = multiplicity ;
    int binNb = -1;
    if (obs >= minvalue) binNb = ((obs - minvalue) / (maxvalue - minvalue)) * nbinstot ;    
    if ( trigger and (binNb >= 0) && (binNb < nbins) ) {
      //std::cout << observable << setw(10) <<  obs << setw(10) << binNb << setw(10)<< npoi  << std::endl;
      binCounts[binNb][1] += npoi;
      binCounts[binNb][2] += 1.;
      binCounts[nbins][1] += npoi;                // <------ used
      binCounts[nbins+1][1] += multiplicity;      // <-------for
      binCounts[nbins+2][1] += 1.;                // <-------normalization
    }
    
    //consider updated bincounts when trigger conditions are satisfied
    if (trigger) {
      subanalysis.incrementNumberOfTriggeredEvents();  //count valid events
      subanalysis.setBinCounts(binCounts);             //put updated binCounts back into subanalysis
    }
  }
};
