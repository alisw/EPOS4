#ifndef PARTICLE_SPECTRA_CENTRALITY_TRIGGER
#define PARTICLE_SPECTRA_CENTRALITY_TRIGGER

#include "Analysis.hpp"

/**
 * @brief transverse impulsion distribution  
 * select the output particles with a given particle id (subparameters[0]) and particle status (0),
 * with a transverse momentum in range[minvalue, maxvalue]  
 * with a rapidity in range[subparameters[1], subparameters[2]]  
 * only if the centrality is in range[subparameters[3],subparameters[4]]
 */

class ParticleSpectraCentralityTrigger : public Analysis {
public:
  using Analysis::Analysis;
  void analyze();
};

#endif
