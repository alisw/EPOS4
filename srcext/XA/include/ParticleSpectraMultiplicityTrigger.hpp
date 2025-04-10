#ifndef PARTICLE_SPECTRA_MULTIPLICITY_TRIGGER
#define PARTICLE_SPECTRA_MULTIPLICITY_TRIGGER

#include "Analysis.hpp"


/**
 * @brief select the output particles with a given particle id (subparameters[0]) 
 * and particle status (0).
 * Other filters are applied according to subanalysis meaning and subparameters values.
 */

class ParticleSpectraMultiplicityTrigger : public Analysis {
public:
  using Analysis::Analysis;
  void analyze();
};

#endif
