#include "Analysis.hpp"
#include "Utils.hpp"

Analysis::Analysis(std::string name, float minvalue, float maxvalue, int nbBins,
                   std::string binning, int normalization, int nbMoreparameters,
                   int nbSubanalyses)
    : name(name), observable(""), minvalue(minvalue), maxvalue(maxvalue),
      numberOfBins(nbBins), binning(binning), normalization(normalization) {
  this->moreparameters.resize(nbMoreparameters);
  this->subanalyses.resize(nbSubanalyses);
}

Analysis::Analysis(std::string name, std::string observable, float minvalue,
                   float maxvalue, int nbBins, std::string binning,
                   int normalization, int nbMoreparameters, int nbSubanalyses)
    : name(name), observable(observable), minvalue(minvalue),
      maxvalue(maxvalue), numberOfBins(nbBins), binning(binning),
      normalization(normalization) {
  this->moreparameters.resize(nbMoreparameters);
  this->subanalyses.resize(nbSubanalyses);
}

Analysis::Analysis(std::string name, std::string observable, float minvalue,
                   float maxvalue, int nbBins, std::string binning,
                   int normalization,
		   std::string histofilename, std::string histoname, 
		   int nbMoreparameters, int nbSubanalyses)
    : name(name), observable(observable), minvalue(minvalue),
      maxvalue(maxvalue), numberOfBins(nbBins), binning(binning),
      normalization(normalization) {
  this->moreparameters.resize(nbMoreparameters);
  this->subanalyses.resize(nbSubanalyses);
  this->initializeImportarray(histofilename, histoname);
}

float Analysis::getMoreparametersValue(unsigned int index) const {
  float value = 0.;
  try {
    value = this->moreparameters.at(index);
  } catch (std::exception &e) {
    std::cerr << std::string(25, '*') << std::endl;
    std::cerr << "Please check the number of more parameters." << std::endl;
    std::cerr << "Size of more parameter array: " << this->moreparameters.size()
              << std::endl;
    std::cerr << "Standard exception: " << e.what() << std::endl;
    std::cerr << std::string(25, '*') << std::endl;
    exit(1);
  }

  return value;
}

float Analysis::getMinvalue() const { return this->minvalue; }

float Analysis::getMaxvalue() const { return this->maxvalue; }

std::string Analysis::getObservable() const { return this->observable; }

unsigned int Analysis::getNumberOfBins() const {
  return (unsigned int)this->numberOfBins;
}

unsigned int Analysis::getNormalization() const {
  return (unsigned int)this->normalization;
}

std::string Analysis::getName() const { return this->name; }

std::vector<std::vector<float>> Analysis::getImportarray() const {
  return this->importarray;
}

void Analysis::setMoreparameters(std::vector<float> moreparametersParam) {
  this->moreparameters = moreparametersParam;
}

std::vector<float> Analysis::getMoreparameters() const {
  return this->moreparameters;
}

unsigned int Analysis::getNumberOfSubanalyses() const {
  return this->subanalyses.size();
}

void Analysis::setSubanalyses(std::vector<Subanalysis> subanalysisParam) {
  this->subanalyses = subanalysisParam;
}

std::vector<Subanalysis> &Analysis::getSubanalyses() {
  return this->subanalyses;
}

void Analysis::addSubparameters(unsigned int subanalysisIndex,
                                std::string meaning,
                                std::vector<float> subparameters) {
  this->subanalyses[subanalysisIndex].setMeaning(meaning);
  this->subanalyses[subanalysisIndex].setSubparameters(subparameters);
}

void Analysis::initializeSubanalysesBinCounts() {
  std::vector<float> binvalues;
  float step, base, value;

  if (this->binning == "lin") {
    step =
        (this->getMaxvalue() - this->getMinvalue()) / this->getNumberOfBins();
    for (unsigned int i = 0; i < this->getNumberOfBins(); i++) {
      value = this->getMinvalue() + (i + 0.5) * step;
      binvalues.push_back(value);
    }
  } else if (this->binning == "log") {
    base = (this->getMaxvalue() / this->getMinvalue());
    for (unsigned int i = 0; i < this->getNumberOfBins(); i++) {
      value = this->getMinvalue() *
              pow(base, ((i + 0.5) / this->getNumberOfBins()));
      binvalues.push_back(value);
    }
  }

  for (Subanalysis &subanalysis : this->subanalyses) {
    subanalysis.createBinCounts(this->getNumberOfBins());
    subanalysis.initializeBinCounts(binvalues);
  }
}

void Analysis::normalize(unsigned int numberOfEvents) {
  unsigned int normalization = this->getNormalization();
  bool normalizeByBinCounts = false;

  // from src/KW/xan.f
  // c.......here normalization.......................................
  // c           see also   "..........fill histogram"
  // c.................................................................
  // c     the norm ( inorm(n) ) is a number hijk which normalizes to:
  // c
  // c  k  0:  * 1
  // c     1:  / number of events
  // c     2:  / number of triggered events
  // c     4:  / bin-counts
  // c
  // c  j  0:  * 1
  // c     1:  / bin-width
  // c
  // c.................................................................
  if (normalization) {

    unsigned int j = normalization / 10;
    unsigned int k = normalization - (j * 10);
    
    for (Subanalysis &subanalysis : this->subanalyses) {
      float factor = 1.;
      switch (j) {
      case 1:
	// binWidth
	factor /=
          (this->getMaxvalue() - this->getMinvalue()) / this->getNumberOfBins();
	break;
      }
      switch (k) {
      case 1:
	factor /= numberOfEvents;
	this->histoweight = numberOfEvents;
	break;
      case 2:
	if (subanalysis.getNumberOfTriggeredEvents() == 0){
	  // set factor to 0 if there is no valid event.
	  factor = 0.;
	}	else {
	  factor /= subanalysis.getNumberOfTriggeredEvents();
	}
	this->histoweight = subanalysis.getNumberOfTriggeredEvents();
	break;
      case 4:
	// binCounts
	normalizeByBinCounts = true;
	this->histoweight = 0.;
	break;
      }
      subanalysis.normalize(factor, normalizeByBinCounts);
    }
  }
  else {
    for (Subanalysis &subanalysis : this->subanalyses) {
      this->histoweight = -1;
    }
  }
}

void Analysis::printBinCounts(unsigned int subanalysisIndex, std::ostream &os) {
  this->subanalyses[subanalysisIndex].printBinCounts(os);
}

unsigned int Analysis::getBinNb(float value) {
  unsigned int binIndex;

  if (this->binning == "lin") {
    binIndex = ((value - this->getMinvalue()) /
                (this->getMaxvalue() - this->getMinvalue())) *
               this->getNumberOfBins();
  } else if (this->binning == "log") {
    binIndex = log(value / this->getMinvalue()) /
               log(this->getMaxvalue() / this->getMinvalue()) *
               this->getNumberOfBins();
  }

  return binIndex;
}

double Analysis::getHistoweight() const { return this->histoweight; }

void Analysis::setHistoweight(double histoweight) {
  this->histoweight = histoweight;
}

void Analysis::initializeImportarray(std::string histofilename, std::string histoname){
  std::string arrayValues = readValuesFromFile(histofilename, histoname);
  this->importarray = readValuesFromString(arrayValues);
}

std::ostream &operator<<(std::ostream &os, Analysis &analysis) {
  os << "Analysis: " << analysis.getName() << std::endl;

  if (analysis.getObservable().length() > 0) {
    os << "Observable: " << analysis.getObservable() << std::endl;
  }

  os << "xmin: " << analysis.getMinvalue()
     << "\txmax: " << analysis.getMaxvalue()
     << "\tnumber of bins: " << analysis.getNumberOfBins()
     << "\tnormalization: " << analysis.getNormalization() << std::endl;
  os << "more parameters : | ";

  for (float parameter : analysis.getMoreparameters()) {
    os << parameter << " | ";
  }
  os << std::endl;

  int n = 0;
  // Analysis::getSubanalyses() returns a Vector of Subanalysis
  for (Subanalysis &subanalysis : analysis.getSubanalyses()) {
    os << "Subanalysis : " << n++ << std::endl;
    os << subanalysis;
  }
  os << std::endl;

  return os;
}
