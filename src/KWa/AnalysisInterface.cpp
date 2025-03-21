#include "Analysis.hpp"
#include "MyAnalysis.hpp"
#include "ParticlePairsCentralityTrigger.hpp"
#include "ParticleSpectraCentralityTrigger.hpp"
#include "ParticleSpectraPercentileTrigger.hpp"
#include "ParticleSpectraMultiplicityTrigger.hpp"
#include "PtSpectra.hpp"
#include "PtSpectraMultTrigger.hpp"
#include "MultiplicityDistributions.hpp"
#include "MultiplicityDependences.hpp"
#include "FlowScalarProductCentralityTrigger.hpp"
#include <cstring>
#include <numeric>
#include <sstream>

std::vector<Analysis *> analysisVector;

extern "C" {
// the following functions are called from Fortran source code

/**
 * @brief create an analysis
 * after parsing the following line with the subroutine aread from src/KW/bas.f
 *  <pre>
 *  beginanalysis
 *   name PtSpectra
 *   minvalue 0.
 *   maxvalue 10.
 *   nrbins   10
 *   binning lin
 *   normalization 1
 *   moreparameters 1
 *     3     !another common parameter
 *   subanalysis 2 !which allows to make similar analyses in parallel
 *   subparameters 3
 *      !First set
 *      120  !particle id
 *      -1.0 !min rapidity
 *       1.0 !max rapidity
 *      !Second set
 *      2130 !particle id
 *      -1.0 !min rapidity
 *       1.0 !max rapidity
 *  endanalysis
 *  </pre>
 *
 * @param name : child class type
 * @param nameLength : analysis name length
 * @param observable : observale name
 * @param observableLength : analysis name length
 * @param minvalue : min X value
 * @param maxvalue : max X value
 * @param nrbins : number of bins
 * @param binning : binning type - linear (lin) or logarithmic (log)
 * @param binningLength : binning type length
 * @param normalization : normalization
 * @param nbOfMoreparameters : number of additional common parameters at least 4
 * @param histofile : histofile name
 * @param histofileLength : length of histofile name
 * @param histogram : histogram name
 * @param histogramLength : length of histogram name
 * @param moreparameters : additional common parameters
 * @param nbOfSubanalyses : number of Subanalysis
 */
void createnamedanalysis_(char *name, int *nameLength, char *observable,
                          int *observableLength, float *minvalue,
                          float *maxvalue, int *nrbins, char *binning,
                          int *binningLength, int *normalization,
			  char *histofile, int *histofileLength,
			  char *histogram, int *histogramLength,
                          int *nbOfMoreparameters, float *moreparameters,
                          int *nbOfSubanalyses) {
  Analysis *analysis = NULL;
  // get analysis name : corresponds to an Analysis child class name (PtSpectra,
  // PtSpectraMultTrigger...)

  std::string analysisName = "";
  for (unsigned int i = 0; i < *nameLength; i++) {
    analysisName = analysisName + name[i];
  }

  std::string observableName = "";
  for (unsigned int i = 0; i < *observableLength; i++) {
    observableName = observableName + observable[i];
  }

  std::string binningType = "";
  for (unsigned int i = 0; i < *binningLength; i++) {
    binningType = binningType + binning[i];
  }

  std::string histofileName = "";
  for (unsigned int i = 0; i < *histofileLength; i++) {
    histofileName = histofileName + histofile[i];
  }

  std::string histogramName = "";
  for (unsigned int i = 0; i < *histogramLength; i++) {
    histogramName = histogramName + histogram[i];
  }

  // each class inherited from the Analysis class must be instanciated here
  // create an Analysis instance from analysis name
  if (analysisName == "PtSpectra") {
    analysis = new PtSpectra(analysisName, *minvalue, *maxvalue, *nrbins, binningType,
                      *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "PtSpectraMultTrigger") {
    analysis = new PtSpectraMultTrigger(analysisName, *minvalue, *maxvalue,
                                        *nrbins, binningType, *normalization,
                                        *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "MyAnalysis") {
    analysis = new MyAnalysis(analysisName, *minvalue, *maxvalue, *nrbins, binning,
                       *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "ParticleSpectraMultiplicityTrigger") {
    analysis = new ParticleSpectraMultiplicityTrigger(
        analysisName, observableName, *minvalue, *maxvalue, *nrbins,
        binningType, *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "ParticleSpectraCentralityTrigger") {
    analysis = new ParticleSpectraCentralityTrigger(
        analysisName, observableName, *minvalue, *maxvalue, *nrbins,
        binningType, *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "ParticleSpectraPercentileTrigger") {
    analysis = new ParticleSpectraPercentileTrigger(
        analysisName, observableName, *minvalue, *maxvalue, *nrbins,
        binningType, *normalization, histofileName, histogramName,
	*nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "MultiplicityDistributions") {
    analysis = new MultiplicityDistributions(
        analysisName, observableName, *minvalue, *maxvalue, *nrbins,
        binningType, *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "ParticlePairsCentralityTrigger") {
    analysis = new ParticlePairsCentralityTrigger(
        analysisName, observableName, *minvalue, *maxvalue, *nrbins,
        binningType, *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "MultiplicityDependences") {
    analysis = new MultiplicityDependences(
        analysisName, observableName, *minvalue, *maxvalue, *nrbins,
        binningType, *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  } else if (analysisName == "FlowScalarProductCentralityTrigger") {
    analysis = new FlowScalarProductCentralityTrigger(
        analysisName, observableName, *minvalue, *maxvalue, *nrbins,
        binningType, *normalization, *nbOfMoreparameters, *nbOfSubanalyses);
  // add a new Analysis instance just before this line
  }

  // if analysis name is a known Analysis child class name
  if (analysis != NULL) {
    // create and fill a c++ vector of common parameter
    std::vector<float> moreparametersVector;
    for (unsigned int i = 0; i < *nbOfMoreparameters; i++) {
      moreparametersVector.push_back(moreparameters[i]);
    }
    // define the subparameters attribute for this Analysis instance
    analysis->setMoreparameters(moreparametersVector);
    // initialize the parameter set binCounts for this analysis
    analysis->initializeSubanalysesBinCounts();
    // add this Analysis instance to the analysisVector : a global variable
    // containing all the Analysys instance
    analysisVector.push_back(analysis);
  } else {
    // the analysis name do not match to a Analysis child class
    std::cout << "unknown analysis type: " << analysisName << std::endl;
    exit(1);
  }
}

/**
 * @brief set common parameters
 *
 * @param nbOfMoreparameters number of common parameters
 * @param moreparameters common parameters array
 */
void setmoreparameters_(int *nbOfMoreparameters, float *moreparameters) {
  // get the last analysis from analysisVector
  Analysis *analysis = analysisVector.back();
  // create a vector of common parameters
  std::vector<float> moreparametersVector;
  for (unsigned int i = 0; i < *nbOfMoreparameters; i++) {
    moreparametersVector.push_back(moreparameters[i]);
  }
  // set this vector as the common parameter attribute
  analysis->setMoreparameters(moreparametersVector);
  // initialize the parameter set binCounts
  analysis->initializeSubanalysesBinCounts();
}

/**
 * @brief set additional parameter to the subanalysis with index
 * subanalysisIndex of the last created analysis
 *
 * @param subanalysisIndex parameter set index
 * @param meaningValue meaning of this additional parameter
 * @param meaningLength length of the string meaningValue
 * @param nbOfSubparameters number of additional parameters
 * @param subparameters array of additional parameters
 */
void addsubparameters_(int *subanalysisIndex, char *meaningValue,
                       int *meaningLength, int *nbOfSubparameters,
                       float *subparameters) {
  std::string meaning = "";
  for (unsigned int i = 0; i < *meaningLength; i++) {
    meaning = meaning + meaningValue[i];
  }

  std::vector<float> subparametersVector;
  for (unsigned int i = 0; i < *nbOfSubparameters; i++) {
    subparametersVector.push_back(subparameters[i]);
  }
  analysisVector.back()->addSubparameters(*subanalysisIndex, meaning,
                                          subparametersVector);
}

/**
 * @brief display all analysis from analysisVector
 * only used for debug
 */
void displayanalysis_() {
  for (Analysis *analysis : analysisVector) {
    std::cout << *analysis << std::endl << std::flush;
  }
}

/**
 * @brief get binCounts values for the rowNb first rows stored in subanalysis
 * with index subanalysisIndex, for analysis with index analysisIndex
 *
 * @param analysisIndex analysis index
 * @param subanalysisIndex parameter set index
 * @param rowNb row number
 * @param values array of rowNb * numberOfValues values
 * @param numberOfValues number of values per row (NR_OF_COLUMNS)
 */
void getbincountsvalues_(int *analysisIndex, int *subanalysisIndex, int *rowNb,
                         float *values, int *numberOfValues) {
  *numberOfValues = NR_OF_COLUMNS;
  try {
    Analysis *analysis = analysisVector.at(*analysisIndex);
    std::vector<Subanalysis> subanalysis = analysis->getSubanalyses();
    std::vector<std::array<float, NR_OF_COLUMNS>> binCounts =
        subanalysis.at(*subanalysisIndex).getBinCounts();

    for (unsigned int i = 0; i < *rowNb; i++) {
      for (unsigned int j = 0; j < *numberOfValues; j++) {
        values[i * (*numberOfValues) + j] = binCounts[i].data()[j];
      }
    }
  } catch (std::exception &e) {
    std::cerr << std::string(25, '*') << std::endl;
    std::cerr << "Please check the number of common parameters." << std::endl;
    std::cerr << "Size of analysis vector: " << analysisVector.size();
    std::cerr << " index of analysis: " << *analysisIndex << std::endl;
    std::cerr << "Size of subanalysis vector: "
              << analysisVector.at(*analysisIndex)->getSubanalyses().size();
    std::cerr << " index of subanalysis: " << *subanalysisIndex << std::endl;
    std::cerr << "Standard exception: " << e.what() << std::endl;
    std::cerr << std::string(25, '*') << std::endl;
    exit(1);
  }
}

/**
 * @brief returns the number of bins for a given analysis
 * the number of bins of a given analysis corresponds to the number of rows of
 * all the arrays stored in the subanalysis of these analysis
 *
 * @param analysisIndex analysis index in analysisVector
 * @return int number of parameter sets
 */
int getnumberofrows_(int *analysisIndex) {
  return analysisVector.at(*analysisIndex)->getNumberOfBins();
}

/**
 * @brief returns the number of sub analyses for a given analysis
 *
 * @param analysisIndex analysis index in analysisVector
 * @return int number of parameter sets
 */
int getnumberofsubanalyses_(int *analysisIndex) {
  return analysisVector.at(*analysisIndex)->getNumberOfSubanalyses();
}

/**
 * @brief get histogram weight for the analysis with analysisIndex in
 * analysisVector
 *
 * @param analysisIndex analysis index in analysisVector
 * @return double histogram weight
 */
double gethistoweight_(int *analysisIndex) {
  return analysisVector.at(*analysisIndex)->getHistoweight();
}

/**
 * @brief delete all Analysis instances from analysisVector
 *
 */
void destroyanalysisvector_() {
  for (Analysis *analysis : analysisVector) {
    delete analysis;
  }
}
}

/**
 * @brief normalize all Analysis instances from analysisVector
 *
 * @param numberOfEvents total number of events
 */
void normalizeAnalysis(const int numberOfEvents) {
  for (unsigned int analysisIndex = 0; analysisIndex < analysisVector.size();
       analysisIndex++) {
    analysisVector.at(analysisIndex)->normalize(numberOfEvents);
  }
}

std::vector<Analysis *> getAnalysesVector() { return analysisVector; }
