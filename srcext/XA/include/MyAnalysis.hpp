#ifndef MYANALYSIS
#define MYANALYSIS

#include "Analysis.hpp"

/**
 * @brief example of custom analysis for transverse impulsion distribution  
 * select the output particles with a given particle id (subparameters[0]) and particle status (0),
 * with a transverse momentum in range[minvalue, maxvalue]  
 * with a rapidity in range[subparameters[1], subparameters[2]]  
 */

class MyAnalysis : public Analysis {
public:
  using Analysis::Analysis;
  void analyze();
};

#endif
