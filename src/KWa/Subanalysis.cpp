#include "Subanalysis.hpp"
#include <cmath>
#include <cstdlib>

Subanalysis::Subanalysis() { this->setNumberOfTriggeredEvents(0); }

Subanalysis::Subanalysis(const int nbOfSubparameters)
    : numberOfTriggeredEvents(0) {
  this->subparameters.resize(nbOfSubparameters);
}

Subanalysis::Subanalysis(const int nbOfBins, const int nbOfSubparameters)
    : numberOfTriggeredEvents(0) {
  this->binCounts.resize(nbOfBins);
  this->subparameters.resize(nbOfSubparameters);
}

Subanalysis::Subanalysis(const std::string meaning, const int nbOfBins,
                         const int nbOfSubparameters)
    : numberOfTriggeredEvents(0), meaning(meaning) {
  this->binCounts.resize(nbOfBins);
  this->subparameters.resize(nbOfSubparameters);
  this->setNumberOfTriggeredEvents(0);
}

Subanalysis::~Subanalysis() {}

float Subanalysis::getSubparametersValue(unsigned int index) const {
  float value = 0.;
  try {
    value = this->subparameters.at(index);
  } catch (std::exception &e) {
    std::cerr << std::string(80, '#') << std::endl;
    std::cerr << "ERROR related to getSubparametersValue." << std::endl;
    std::cerr << "Size of subparameters array: " << this->subparameters.size() << std::endl;
    std::cerr << "Maximal possible index: " <<  this->subparameters.size() -1 << std::endl;
    std::cerr << "Current index: " << index  << std::endl;
    //std::cerr << "Standard exception: " << e.what() << std::endl;
    std::cerr << std::string(80, '#') << std::endl;
    exit(1);
  }

  return value;
}

void Subanalysis::setSubparameters(
    const std::vector<float> subparametersParam) {
  // perform a deep copy
  this->subparameters.assign(subparametersParam.begin(),
                             subparametersParam.end());
}

std::vector<float> Subanalysis::getSubparameters() const {
  return this->subparameters;
}

std::vector<std::array<float, NR_OF_COLUMNS>>
Subanalysis::getBinCounts() const {
  return this->binCounts;
}

void Subanalysis::setBinCounts(
    std::vector<std::array<float, NR_OF_COLUMNS>> &binCounts) {
  this->binCounts = binCounts;
}

std::string Subanalysis::getMeaning() const { return this->meaning; }

void Subanalysis::setMeaning(std::string meaning) { this->meaning = meaning; }

void Subanalysis::createBinCounts(unsigned int nbOfBins) {
  this->binCounts.resize(nbOfBins);
}

void Subanalysis::initializeBinCounts(std::vector<float> binvalues) {
  int i = 0;
  for (float binvalue : binvalues) {
    this->binCounts[i++] = {binvalue, 0., 0.};
  };
}

void Subanalysis::update(unsigned int binNb, float value, float error) {
  this->binCounts[binNb][1] += value;
  this->binCounts[binNb][2] += error;
}

void Subanalysis::reset() {
  for (std::array<float, NR_OF_COLUMNS> &row : this->binCounts) {
    row[1] = 0.;
    row[2] = 0.;
  }
}

void Subanalysis::printBinCounts(std::ostream &os) {
  os << "binCounts " << NR_OF_COLUMNS << std::endl;
  os << std::scientific;
  for (std::array<float, NR_OF_COLUMNS> row : this->getBinCounts()) {
    for (float value : row) {
      os << "\t" << value;
    }
    os << std::endl;
  }
  os << "end binCounts " << std::endl;
}

float *Subanalysis::getBinCountsValues(int rowNumber) {
  return this->getBinCounts()[rowNumber].data();
}

void Subanalysis::setNumberOfTriggeredEvents(
    unsigned int numberOfTriggeredEvents) {
  this->numberOfTriggeredEvents = numberOfTriggeredEvents;
}

unsigned int Subanalysis::getNumberOfTriggeredEvents() {
  return this->numberOfTriggeredEvents;
}

void Subanalysis::incrementNumberOfTriggeredEvents() {
  this->numberOfTriggeredEvents += 1;
}

void Subanalysis::normalize(float factor, bool normalizeByBinCounts) {
  float row_tmp;
  for (std::array<float, NR_OF_COLUMNS> &row : this->binCounts) {
    if (normalizeByBinCounts and row[2]) {
      row[1] /= row[2];
    } else {
      row[1] *= factor;
      if (row[2]) {
        row[2] = row[1] / sqrt(row[2]);
      }
    }
  }
}

std::ostream &operator<<(std::ostream &os, const Subanalysis &subanalysis) {
  os << "sub parameters : | ";
  for (float parameter : subanalysis.getSubparameters()) {
    os << parameter << " | ";
  }
  os << std::endl;

  os << "BinCounts" << std::endl;
  for (std::array<float, NR_OF_COLUMNS> row : subanalysis.getBinCounts()) {
    for (float value : row) {
      os << value << "\t| ";
    }
    os << std::endl;
  }
  os << std::endl;

  return os;
}
