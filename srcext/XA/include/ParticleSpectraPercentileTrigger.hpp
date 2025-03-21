#ifndef PARTICLE_SPECTRA_PERCENTILE_TRIGGER
#define PARTICLE_SPECTRA_PERCENTILE_TRIGGER

#include "Analysis.hpp"


/**
 * @brief transverse impulsion distribution  
 * select the output particles with a given particle id (subparameters[0]) and particle status (0),
 * with a transverse momentum in range[minvalue, maxvalue]  
 * with a rapidity in range[subparameters[1], subparameters[2]]  
 * only if the multiplicity is in range[subParameters[3],subParameters[4]]
 */

class ParticleSpectraPercentileTrigger : public Analysis {
public:
  using Analysis::Analysis;
  void analyze();
};

#endif
