//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License
//  version 3 or later (See COPYING file for the text of the licence)
//

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Units.h"
#include "HepMC3/Version.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>


#if defined HEPMC3_USE_COMPRESSION && defined HEPMC3_Z_SUPPORT
#include "HepMC3/WriterGZ.h"
HepMC3::WriterGZ<HepMC3::WriterAscii> *writer;
#else
HepMC3::WriterAscii *writer;
#endif

std::shared_ptr<HepMC3::GenRunInfo> run;

extern "C" void openhepmc_(const char *filename);

extern "C" void fillhepmc_(int *_iextree, int *_ihepmc3, int *_nevt, float *_eng, float *_dyframe,
                    int *_iprojZ, int *_iprojA, int *_itargZ, int *_itargA,
                    int *_nhard, int *_ncoll, int *_npartproj, int *_nparttarg,
                    int *_nspecprojp, int *_nspecprojn, int *_nspectargp, int *_nspectargn, 
                    int *_np, float *_bim, float *_sigtot,
                    int *_id, int *_ist, int *_ity, int *_ior, float *_px, float *_py,
                    float *_pz, float *_m, float *_x, float *_y, float *_z, float *_t,
                    int *_record_mode, float *_tau_decay, int *_record_id_nb,
                    int *_record_id_list);

extern "C" void closehepmc_();

extern "C" void id_epostopdg_(int *_idepos, int *_idpdg);

extern "C" void idtau_(int *_id, float *_p4, float *_p5, float *_tau);


/**
 * @fn void openhepmc_(char* filename)
 * @brief Function creating the HepMC event object
 *
 * @param filename HepMC filename
 */

void openhepmc_(const char *filename) {
  run = std::make_shared<HepMC3::GenRunInfo>();
  run->set_weight_names({"first_weight"});

  // generator information
  std::string version;
  std::string generator_version("4.0.0");
  std::string version_filename = std::string(std::getenv("EPO")) + std::string("VERSION.txt");
  std::ifstream version_file(version_filename);
  if (version_file.is_open()) {
    while (getline(version_file, version)) {
      if (version.length() > 0) {
	generator_version = version;
      };
    };
    version_file.close();
  } else {
    std::cout << "Unable to open file ";
    std::cout << version_filename << std::endl;
    std::cout << "set EPOS default version to 4.0.0" << std::endl;
  }
  
  std::string generator_name("EPOS");
  struct HepMC3::GenRunInfo::ToolInfo generator={generator_name, 
						 generator_version,
						 std::string("used generator")};
  run->tools().push_back(generator);

  // configuration file information
  std::string hepmc_filename(filename);
  std::string prefix("z-");
  std::string suffix(".hepmc");
  std::string option_filename = hepmc_filename.substr(hepmc_filename.find_last_of("/") + prefix.length() + 1,std::string::npos);
  std::size_t found_suffix = option_filename.rfind(suffix);
  if (found_suffix != std::string::npos)
    option_filename.replace(found_suffix, suffix.length(), ".optns");
  struct HepMC3::GenRunInfo::ToolInfo config={option_filename,"1.0",std::string("configuration file")};
  run->tools().push_back(config);

#if defined HEPMC3_USE_COMPRESSION && defined HEPMC3_Z_SUPPORT
  std::string hepmc_gz_filename = hepmc_filename + ".gz";
  writer = new HepMC3::WriterGZ<HepMC3::WriterAscii>(hepmc_gz_filename, run);
#else
  writer = new HepMC3::WriterAscii(filename, run);
#endif

  // writer = new HepMC3::WriterAscii(filename);
  if (!writer) {
    std::cout << " FAILED!!!! " << std::endl;
  }
  std::cout << "*** HepMC 3 ***" << std::endl;
  std::cout << "*** " << filename << " open***" << std::endl;
}

/**
 * @fn double get_mass(int A, int Z)
 * @brief Function assigning the mass to the colliding ions
 *
 * @param A mass number
 * @param Z atomic number
 */
double get_mass(int A, int Z) {
  double m = 0;
  //  U :
  if (A == 238 && Z == 92)
    m = 221.742894;
  //  Pb :
  if (A == 208 && Z == 82)
    m = 193.729016;
  //  Au :
  if (A == 197 && Z == 79)
    m = 183.473190;
  //  Xe :
  if (A == 129 && Z == 54)
    m = 120.074038;
  //  Cu :
  if (A == 63 && Z == 29)
    m = 58.6185461;
  //  d :
  if (A == 2 && Z == 1)
    m = 1.87612384;
  //  p  :
  if (A == 1 && Z == 1)
    m = 0.93827208;
  return m;
}

/**
 * @fn void fillhepmc_(int *_iextree, int *_ihepmc3, int *_nevt, float *_eng, 
 float *_dyframe, int *_iprojZ, int *_iprojA, int *_itargZ, int *_itargA, int 
 *_nhard, int *_ncoll, int *_npartproj, int *_nparttarg, int *_nspecprojp, int 
 *_nspecprojn, int *_nspectargp, int *_nspectargn, int *_np, float *_bim, float 
 *_sigtot, int *_id, int *_ist, int *_ity, int *_ior, float *_px, float *_py, 
 float *_pz, float *_e, float *_x, float *_y, float *_z, float *_t, int 
 *_record_mode, float *_tau_decay, int* _record_id_nb, int* _record_id_list)

 * @brief Function filling the HepMC event
 *
 * @param _iextree root tree structure flag : if you want to use the new tree
 structure, set _iextree to 1, in order to use the old tree structure, _iextree
 has to be set to 0
 * @param _ihepmc3 root tree/hepmc structure flag : if you want to use the new 
 root tree and hepmc3 event structures, set _ihepmc3 to 1; in order to use the 
 old tree and hepmc2 event strucures, _ihepmc3 has to be set to 0
 * @param _nevt event number
 * @param _eng center-of-mass energy in GeV
 * @param _dyframe center-of-mass rapidity compared to the lab frame
 * @param _iprojZ atomic number projectile
 * @param _iprojA mass number projectile
 * @param _itargZ atomic number target
 * @param _itargA mass number target
 * @param _nhard number of elementary hard parton-parton scatterings (recorded
 only if _iextree == 1)
 * @param _ncoll number of binary collisions according to Glauber (recorded only
 if _iextree == 1)
 * @param _npartproj number of projectile's nucleons participants according to
 Glauber (recorded only if _iextree == 1)
 * @param _nparttarg number of target's nucleons participants according to
 Glauber (recorded only if _iextree == 1)
 * @param _nspecprojp number of projectile's spectators protons according to EPOS 
 (recorded only if _iextree == 1)
 * @param _nspecprojn number of projectile's spectators neutrons according to EPOS 
 (recorded only if _iextree == 1)
 * @param _nspectargp number of target's spectators protons according to EPOS 
 (recorded only if _iextree == 1)
 * @param _nspectargn number of target's spectators neutrons according to EPOS 
 (recorded only if _iextree == 1)
 * @param _np number of particles in the event
 * @param _bim impact parameter
 * @param _sigtot corresponding pp cross-section (recorded only if _iextree ==
 1)
 * @param _id particle id : see tables in file KWt/idt.dt : first column for
 EPOS IDs, second one for PDG IDs
 * @param _ist particle status (hadron last generation (0) or not (1); other
 numbers refer to partons, Pomerons, etc)
 * @param _ity type of particle origin (20-29 from soft strings, 30-39 from hard
 strings, 40-59 from remnants, 60 from fluid)
 * @param _ior index of father
 * @param _px p_x of particle
 * @param _py p_y of particle
 * @param _pz p_z of particle
 * @param _m mass of particle
 * @param _x x component of formation point
 * @param _y y component of formation point
 * @param _z z component of formation point
 * @param _t formation time
 * @param _record_mode defines which particle selection is applied (1 : only FS
 particles || 2 : FS particles + parents & vertex for weak & EM decayed
 particles)
 * @param _tau_decay limit lifetime to distinguish strong from weak/EM decay
 * @param _record_id_nb number of particles ID with tau < tau_decay but for
 which decay record is needed
 * @param _record_id_list list of particles ID with tau < tau_decay but for
 which decay record is needed

 */
void fillhepmc_(int *_iextree, int *_ihepmc3, int *_nevt, float *_eng, float *_dyframe,
                int *_iprojZ, int *_iprojA, int *_itargZ, int *_itargA,
                int *_nhard, int *_ncoll, int *_npartproj, int *_nparttarg,
		int *_nspecprojp, int *_nspecprojn, int *_nspectargp, int *_nspectargn, 
                int *_np, float *_bim, 
                float *_sigtot, int *_id, int *_ist, int *_ity, int *_ior,
                float *_px, float *_py, float *_pz, float *_m, float *_x,
                float *_y, float *_z, float *_t, int *_record_mode,
                float *_tau_decay, int *_record_id_nb, int *_record_id_list) {
  //  ***********************************  //
  //  ****** EVENT INITIALISATION *******  //
  //  ***********************************  //

  //  *** Event definition ***  //
  HepMC3::GenEvent evt(run, HepMC3::Units::GEV, HepMC3::Units::MM); //  create event
  evt.set_event_number(*_nevt);                                //  event number

  //  *** Cross-section ***  //
  //  Generate the cross-section object
  HepMC3::GenCrossSectionPtr cross_section =
      std::make_shared<HepMC3::GenCrossSection>();
  std::vector<double> xs, xs_err;
  xs_err = {0};
  if (*_iextree == 1) {    //  provide value only if new tree used (from EPOS 4)
    xs = {*_sigtot * 1e9}; //  mb (EPOS) -> pb (HepMC)
  } else {
    xs = {0.0};
  }
  cross_section->set_cross_section(xs, xs_err);
  evt.add_attribute("GenCrossSection", cross_section);

  // *** Defining the beam particles *** //
  int start_particle_loop = 0;

  //  * A,Z *  //
  int A_proj = *_iprojA;
  int Z_proj = *_iprojZ;
  int A_targ = *_itargA;
  int Z_targ = *_itargZ;

  //  * Projectile *  //
  int id_proj = 99;
  double mass_beam1, energy_beam1;

  //  Ion case
  if (A_proj > 1) {
    //  ID
    id_proj = 1000000000 + Z_proj * 10000 + A_proj * 10;
    //  Mass
    mass_beam1 = get_mass(A_proj, Z_proj);
    //  Energy
    energy_beam1 = *_eng * 0.5 * A_proj;
    //  Other cases (single hadrons/leptons...)
  } else {
    //  ID
    id_epostopdg_(&_id[0], &id_proj);
    //  Mass
    mass_beam1 = _m[0];
    //  Energy
    energy_beam1 = *_eng * 0.5;
    //  Starting particle loop index later, to avoid double-counting
    start_particle_loop = 1;
  }

  //  * Target *  //
  int id_targ = 99;
  double mass_beam2, energy_beam2;

  //  Ion case
  if (A_targ > 1) {
    //  ID
    id_targ = 1000000000 + Z_targ * 10000 + A_targ * 10;
    //  Mass
    mass_beam2 = get_mass(A_targ, Z_targ);
    //  Energy
    energy_beam2 = *_eng * 0.5 * A_targ;
    //  Other cases (single hadrons/leptons...)
  } else {
    //  ID
    id_epostopdg_(&_id[1], &id_targ);
    //  Mass
    mass_beam2 = _m[1];
    //  Energy
    energy_beam2 = *_eng * 0.5;
    //  Starting particle loop index later, to avoid double-counting
    start_particle_loop = 2;
  }

  //  * Boost in different frame if necessary *  //
  if (*_dyframe != 0.0) {
    if (*_nevt == 1) {
      std::cout << "  Particles shifted from CoM to lab with dy=" << -*_dyframe
                << std::endl;
    }
    energy_beam1 = cosh(*_dyframe) * energy_beam1 -
                   sinh(*_dyframe) * sqrt(energy_beam1 * energy_beam1 -
                                          mass_beam1 * mass_beam1);
    energy_beam2 = cosh(*_dyframe) * energy_beam2 -  
                   sinh(*_dyframe) * (-1) 
                   *sqrt(energy_beam2 * energy_beam2 - mass_beam2 * mass_beam2);
  }

  //  * Create the corresponding particle objects *  //
  HepMC3::GenParticlePtr gp_beam1 = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(
          0, 0, sqrt(energy_beam1 * energy_beam1 - mass_beam1 * mass_beam1),
          energy_beam1),
      id_proj, 4);
  HepMC3::GenParticlePtr gp_beam2 = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(
          0, 0,
          (-1) * sqrt(energy_beam2 * energy_beam2 - mass_beam2 * mass_beam2),
          energy_beam2),
      id_targ, 4);

  //  * Assign the masses *  //
  gp_beam1->set_generated_mass(mass_beam1);
  gp_beam2->set_generated_mass(mass_beam2);

  //  * Attach them to the event object *  //
  evt.set_beam_particles(gp_beam1, gp_beam2);

  //  *** Zero Vertex ***  //
  HepMC3::GenVertexPtr v0 =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(0, 0, 0, 0));
  evt.add_vertex(v0);

  //  * Adding colliding nuclei to vertex 0 *  //
  v0->add_particle_in(gp_beam1);
  v0->add_particle_in(gp_beam2);

  //  * Set barcodes for beams *  //
  // HepMC2 evt.set_signal_process_id(v0->id());
  // HepMC3 evt.add_attribute("signal_process_id",
  //		    std::make_shared<HepMC3::IntAttribute>(v0->id()));
  // HepMC2 evt.set_signal_process_vertex(v0);
  // HepMC3 evt.add_attribute("signal_vertex_id",
  //		    std::make_shared<HepMC3::IntAttribute>(v0->id()));

  //  * Check interaction boson *  //
  if (_id[2] == 10 and _ist[2] == 1) {
    // Starting particle loop index later, to ignore virtual photon in lepton
    // collisions
    start_particle_loop = 3;
  }

  //  *** Heavy ion information ***  //
  if (A_proj > 1 || A_targ > 1) {
    HepMC3::GenHeavyIonPtr hi = std::make_shared<HepMC3::GenHeavyIon>();
    hi->impact_parameter = *_bim; //  impact parameter
    //  Provide following info only if new tree used (from EPOS 4)
    if (*_iextree == 1) {
      hi->Ncoll_hard = *_nhard;         //  number of hard interactions
      hi->Ncoll = *_ncoll;              //  number of Glauber nucleon collisions
      hi->Npart_proj = *_npartproj;     //  number of participants in projectile
      hi->Npart_targ = *_nparttarg;     //  number of participants in target
      hi->spectator_protons = *_nspecprojp + *_nspectargp;  //  number of protons spectators
      hi->spectator_neutrons = *_nspecprojn + *_nspectargn; //  number of neutrons spectators
      if (*_ihepmc3 == 1) {
        hi->Nspec_proj_p = *_nspecprojp; //  number of projectile's protons spectators
        hi->Nspec_proj_n = *_nspecprojn; //  number of projectile's neutrons spectators
        hi->Nspec_targ_p = *_nspectargp; //  number of target's protons spectators
        hi->Nspec_targ_n = *_nspectargn; //  number of target's neutrons spectators
        hi->sigma_inel_NN = *_sigtot; //  inelastic NN cross-section (in mb)
      }
    }
    evt.set_heavy_ion(hi);
  }

  //  ********************************  //
  //  *** DEFINE GENERAL VARIABLES ***  //
  //  ********************************  //
  //  Total number of particles
  int np = *_np;

  //  PDG ID variable
  int idpdg;

  //  Particle lifetime
  float tau;

  //  Limit lifetime to distinguish strong from weak/EM decay
  float tau_decay = *_tau_decay;
  
  //  EPOS status variables for FS and unstable particles respectively
  int ist1 = 0;
  int ist2 = 1;

  //  Define the "used" list initialised with -1 :
  //  - "-1" = never used
  //  -  "0" = used as vtx output but NOT final-state
  //  -  "1" = final-state OR used as vtx output + input
  int used[np];
  for (int i = 0; i < np; i++)
    used[i] = -1;

  //  Define a vector of pointers which will contain decay children that will
  //  decay again
  // (so avoid declaring twice the same particle)
  std::vector<HepMC3::GenParticlePtr> temp_ptl;

  //  Default energy and mass to have gamma = 1 in idtau_(..)
  float p4 = 1.0, p5 = 1.0;

  //  ***** Selection of the EPOS particles written in the HepMC event *****  //
  switch (*_record_mode) {

  case 0:
    // ****************************************************************************************** //
    // ***** All particles of the event (except mediator bosons and pomerons, strings, etc) ***** //
    // ****************************************************************************************** //
    
    if(*_record_mode == 0) {
      if (*_nevt == 1) {
        // Printing out record mode
        std::cout << "  HepMC record mode: FULL" << std::endl;
      }
    
      // Set limit lifetime to record particles to 0, to record all of them
      tau_decay = 0;
    }

  case 1:
    //  ************************************************  //
    //  ***** Final-state particles only (ist = 0) *****  //
    //  ************************************************  //
   
    if(*_record_mode == 1) {
      if (*_nevt == 1) {
        // Printing out record mode
        std::cout << "  HepMC record mode: FINAL_STATE" << std::endl;
      }

      // Set limit lifetime to record particles to 1s, to record only final-state
      tau_decay = 1;
    }
    
  case 4:
    //  *****************************************************************************************************  //
    //  ***** HACAS FULL: Final-state particles + decays (weak, EM & exceptions) if there were no hacas *****  //
    //  *****************************************************************************************************  //
   
    if(*_record_mode == 4) {
      if (*_nevt == 1) {
        // Printing out record mode
        std::cout << "  HepMC record mode: WITHOUT_HACAS" << std::endl;
      }

      // Modify ist1 & ist2 to consider particles from alternative scenario without hadronic cascades
      ist1 = 8; // status code for FS particle in alternative scenario 
      ist2 = 6; // status code for unstable particle in alternative scenario 
    }
    
  case 2:
    // ****************************************************************** //
    // ***** Final-state particles + decays (weak, EM & exceptions) ***** //
    // ****************************************************************** //

    if (*_nevt == 1) {
      if(*_record_mode == 2) {
        // Printing out record mode
        std::cout << "  HepMC record mode: DECAYS" << std::endl;
      }
    
      // Printing out limit lifetime to record particles
      if(*_record_mode >= 2){
        // Printing out record mode
        std::cout << "    -> tau_decays = " << tau_decay << " s" << std::endl;
      }

      // Printing out IDs of exceptions which are recorded whatever the conditions are
      if(*_record_id_nb > 0) {
          std::cout << "    -> IDs of particles recorded as exceptions: " ;
          // Looping over exceptions IDs 
          for (int e = 0; e < *_record_id_nb; e++) {
            std::cout << *(_record_id_list + e);
            if(e < *_record_id_nb - 1) {
              std::cout << " / ";
            }
          }
          std::cout << std::endl;
      }
    }

    // *** Loop over the EPOS particle list *** //
    for (int p = start_particle_loop; p < np; p++) {

      // ***** For weak or EM decayed particles ***** //
      if (_ist[p] == ist2 && used[p] != 1) {

        //  *** Determines the lifetime*** //
        idtau_(&_id[p], &p4, &p5, &tau);
        //  Calculate tau in seconds from c*tau in fm
        tau = tau * 1.e-15 / 2.99792458e8;

        //  *** Determines the decay type from parent's lifetime ***  //
        //  To mark if decay is strong or weak/EM
        bool weak_EM = false;
        //  Evaluate the decay type (EXCEPT FOR INTERACTING NUCLEONS)
        if (tau > tau_decay && p > A_targ + A_proj - 1)
          weak_EM = true;

        //  *** Record the decay if particle in exceptions list ***  //
        for (int e = 0; e < *_record_id_nb; e++) {
          if (abs(_id[p]) == *(_record_id_list + e)) {
            weak_EM = true;
          }
        }

        // *** STRONG DECAY + INTERACTING NUCLEON + used[p] != 0 -> ignored ***
        // //
        if (weak_EM == false && used[p] == -1) {

          //  Mark parent as 'read'
          used[p] = 1;
          continue;

          //  *** WEAK / EM DECAY + EXCEPTIONS -> recorded ***  //
        } else {
          //  *** Parent's information ***  //
          //  Translate ID
          idpdg = 99;
          id_epostopdg_(&_id[p], &idpdg);
          // Status (2 for 'decayed particle', by default)
          int hep_status = 2;
          // Assignate entry variables (for precision concerns, issues if entry
          // variables used directly)
          double m = _m[p];
          double px = _px[p];
          double py = _py[p];
          double pz = _pz[p];
	  
          // Calculate E
          double mt = sqrt(m * m + px * px + py * py);
          double E = mt * cosh(asinh(pz / mt));
          //  Shift E and pz in the lab frame, if necessary
          if (*_dyframe != 0.0) {
            double E_lab = E * cosh(*_dyframe) - pz * sinh(*_dyframe);
            double pz_lab = pz * cosh(*_dyframe) - E * sinh(*_dyframe);
            E = E_lab;
            pz = pz_lab;
          }
          //  K0(b) (EXCEPTION): only record directly its K0s/L child
          int p_K0 =
              p; //  record the index in order to come back to it at the end
          bool K0 =
              false; //  to indicate if child is a K0 that has been ignored
          if (abs(idpdg) == 311) {
            // Indicate that a K0(b) has been found
            K0 = true;
            // Search for the child K0s/L
            for (int K0sl = p + 1; K0sl < np; K0sl++) {
              if (_ior[K0sl] == p + 1) {
                // Mark the K0(b) as 'used'
                used[p] = 1;
                // Change index to the K0s/L one (for particle information)
                p = K0sl;
                // Translate ID of the K0s/L
                id_epostopdg_(&_id[K0sl], &idpdg);
                // Check status of the K0s/L (change to 1 if FS) 
                if(_ist[p] == ist1) {
                  hep_status = 1;
                }
                // Break the loop
                break;
              }
            }
          }
          // Define the particle object and its mass
          HepMC3::GenParticlePtr parent = std::make_shared<HepMC3::GenParticle>(
              HepMC3::FourVector(px, py, pz, E), idpdg, hep_status);            

          // Define the mass of the particle
          parent->set_generated_mass(m);

          // *** If parent is NOT originating from a decay *** //
          if (used[p] == -1) {
            // Attach the parent as a child of the original vertex
            v0->add_particle_out(parent);
          // *** If parent is originating from a decay *** //
          } else {
            // Search for it in the temporary storing list 'temp_ptl'
            for (std::size_t t = 0; t < temp_ptl.size(); t++) {
              // Test if parent is the particle checked in 'temp_ptl' 
              // (compare mass, energy and px)
              if (m == temp_ptl[t]->generated_mass() &&
                  E == temp_ptl[t]->momentum().e() &&
                  px == temp_ptl[t]->momentum().px()) {
                // Test if ID is the same
                if (idpdg != temp_ptl[t]->pdg_id()) {
                  std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°"
                            << evt.event_number() << ") :" << std::endl;
                  std::cout << "     Found correspondence between current "
                               "particle and a decay child (PDG_ID = "
                            << temp_ptl[t]->pdg_id()
                            << "), but ID doesn't match" << std::endl;
                  std::cout << "\t Current particle : EPOS_ID = " << _id[p]
                            << "  /  corresponding PDG_ID = " << idpdg
                            << " (from id_epostopdg_)" << std::endl;
                  std::cout
                      << "     -> Recorded as the same particle anyway, but "
                         "CHECK id_epostopdg_ (KW/ids.f) for this ID\n"
                      << std::endl;
                }
                // Assign 'parent' to its already created corresponding object
                parent = temp_ptl[t];
                // Remove the object from the 'temp_ptl' list
                temp_ptl.erase(temp_ptl.begin() + t);
              }
            }
          }
          // K0(b) (EXCEPTION)
          if (K0 == true) {
            // If K0s/L is FS particle...
            if(hep_status == 1) {
              // ...mark it as 'FS / used for output + input'
              used[p] = 1;
              // ...re-assign 'p' counter to the originally found K0(b) index to 
              // resume the particle loop at initial position 
              p = p_K0;
              // ...skip to next particle to avoid useless search for children 
              continue;
            }
          }

          // *** Search for the children *** //
          // Marker used to indicate if vertex already defined
          bool vertex = false;
          // Marker used to indicate when to stop reading
          bool family = false;
          // Create the pointer for secondary vertex
	  auto v = std::make_shared<HepMC3::GenVertex>();
          // Loop over the other particles
          for (int c = p + 1; c < np; c++) {
            // Check if child from particle p
            if (_ior[c] == p + 1) {
              // Indicate that a child is found
              family = true;
              // *** Secondary vertex *** //
              if (vertex == false) {
                // Convert units from fm (EPOS) to mm (HepMC)
                double x = _x[c] * 1.0e-12;
                double y = _y[c] * 1.0e-12;
                double z = _z[c] * 1.0e-12;
                double t = _t[c] * 1.0e-12;
                //  Define the vertex
                v->set_position(HepMC3::FourVector(x, y, z, t));
                //  Attach it to the event
                evt.add_vertex(v);
                //  Mark vertex as defined
                vertex = true;
                //  *** Attach the parent to the secondary vertex ***  //
                v->add_particle_in(parent);
              }

              //  *** Child's information ***  //
              //  Translate ID
              idpdg = 99;
              id_epostopdg_(&_id[c], &idpdg);
              //  Assign entry variables (for precision concerns, issues if
              //  entry variables used directly)
              double m = _m[c];
              double px = _px[c];
              double py = _py[c];
              double pz = _pz[c];
	      
              //  Calculate E
              double mt = sqrt(m * m + px * px + py * py);
              double E = mt * cosh(asinh(pz / mt));
              //  Shift E and pz in the lab frame, if necessary
              if (*_dyframe != 0.0) {
                double E_lab = E * cosh(*_dyframe) - pz * sinh(*_dyframe);
                double pz_lab = pz * cosh(*_dyframe) - E * sinh(*_dyframe);
                E = E_lab;
                pz = pz_lab;
              }
              // Define the particle object for child
              HepMC3::GenParticlePtr child;
              //  K0(b) (EXCEPTION): only record directly its K0s/L child
              int c_K0 =
                  c; //  record the index in order to come back to it at the end
              bool K0 =
                  false; //  to indicate if child is a K0 that has been ignored
              if (abs(idpdg) == 311) {
                // Indicate that a K0(b) has been found
                K0 = true;
                // Search for the child K0s/L
                for (int K0sl = c + 1; K0sl < np; K0sl++) {
                  if (_ior[K0sl] == c + 1) {
                    // Mark the K0(b) as 'used'
                    used[c] = 1;
                    // Change index to the K0s/L one (for particle information)
                    c = K0sl;
                    // Translate ID of the K0s/L
                    id_epostopdg_(&_id[K0sl], &idpdg);
                    // Break the loop
                    break;
                  }
                }
              }

              if (_ist[c] == ist1) {
                // If FS particle : set hep_status = 1
                child = std::make_shared<HepMC3::GenParticle>(
                    HepMC3::FourVector(px, py, pz, E), idpdg, 1);
                // Mark child as 'final-state'
                used[c] = 1;
              } else {
                // If not FS particle (to decay later) : keep the object until
                // used as parent + set hep_status = 2
                temp_ptl.push_back(std::make_shared<HepMC3::GenParticle>(
                    HepMC3::FourVector(px, py, pz, E), idpdg, 2));
                // Assign 'child' to its corresponding stored object to be
                // re-used
                child = temp_ptl[temp_ptl.size() - 1];
                // Mark child as 'used as input but NOT final-state'
                used[c] = 0;
              }
              // Define the mass of the child
              child->set_generated_mass(_m[c]);
              // *** Attach the child to the secondary vertex *** //
              v->add_particle_out(child);

              // K0(b) (EXCEPTION)
              if (K0 == true) {
                // Re-assign 'c' counter to the originally found K0(b) index to 
                // resume the particle loop at initial position 
                c = c_K0;
              }

              // Signal if ID correspondance not found with 'id_epostopdg_'
              if (idpdg == 99) {
                std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°"
                          << evt.event_number() << " / vertex n°"
                          << child->production_vertex()->id()
                          << ") :" << std::endl;
                std::cout << "\t Cannot find correspondent PDG_ID for decay "
                             "child particle with EPOS_ID = "
                          << _id[c] << " -> Recorded with PDG_ID = 99"
                          << std::endl;
                std::cout
                    << "     -> CHECK id_epostopdg_ (KW/ids.f) for this ID\n"
                    << std::endl;
              }

              // Check if children have been found already
            } else if (_ior[c] != p + 1 && family == true) {
              // Break the loop if children already found
              // (uses the fact that children from same parent always follow
              // each other in the list)
              break;

              // Read the next particle if not a child
            } else {
              continue;
            }
          }

          // *** No children found : ERROR MESSAGE *** //
          if (family == false) {
            std::cout
                << "\n /!\\ NO CHILDREN FOUND DESPITE ist=1 FOR PARTICLE n°"
                << p << " (ID = " << _id[p] << ") \n"
                << std::endl;
            // Set hep_status = 1 then (FS particle)
            parent->set_status(12);
          }

          // Mark parent as 'FS / used for output + input'
          used[p] = 1;

          // K0(b) (EXCEPTION)
          if (K0 == true) {
            // Re-assign 'p' counter to the originally found K0(b) index to 
            // resume the particle loop at initial position 
            p = p_K0;
          }
        }
      // ***** For other FS particles (not from weak or EM decay) ***** //
      } else if (_ist[p] == ist1 && used[p] == -1) {
        // Translate ID
        idpdg = 99;
        id_epostopdg_(&_id[p], &idpdg);
        // Assign entry variables (for precision concerns, issues if entry
        // variables used directly)
        double m = _m[p];
        double px = _px[p];
        double py = _py[p];
        double pz = _pz[p];
	
        // Calculate E
        double mt = sqrt(m * m + px * px + py * py);
        double E = mt * cosh(asinh(pz / mt));
        // Shift E and pz in the lab frame, if necessary
        if (*_dyframe != 0.0) {
          double E_lab = E * cosh(*_dyframe) - pz * sinh(*_dyframe);
          double pz_lab = pz * cosh(*_dyframe) - E * sinh(*_dyframe);
          E = E_lab;
          pz = pz_lab;
        }
        // Define the particle object and its mass
        HepMC3::GenParticlePtr child = std::make_shared<HepMC3::GenParticle>(
            HepMC3::FourVector(px, py, pz, E), idpdg, 1);
        // Define the mass of the particle
        child->set_generated_mass(m);
        // Mark it as 'read'
        used[p] = 1;
        // *** Attach the particle to the original vertex *** //
        v0->add_particle_out(child);

        // Signal if ID correspondance not found with 'id_epostopdg_'
        if (idpdg == 99) {
          std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°"
                    << evt.event_number() << " / vertex n°"
                    << child->production_vertex()->id() << ") :" << std::endl;
          std::cout << "\t Cannot find correspondent PDG_ID for decay child "
                       "particle with EPOS_ID = "
                    << _id[p] << " -> Recorded with PDG_ID = 99" << std::endl;
          std::cout << "     -> CHECK id_epostopdg_ (KW/ids.f) for this ID\n"
                    << std::endl;
        }
      }
    }
    break;

  case 3:
    // ******************************************** //
    // ***** Hadrons before hadronic cascades ***** //
    // ******************************************** //

    // Printing out record mode
    std::cout << "  HepMC record mode: BEFORE_HACAS" << std::endl;

    // *** Loop over the EPOS particle list *** //
    for (int p = start_particle_loop; p < np; p++) {
      if (_ist[p] == 3 && used[p] == -1) {
        // Translate ID
        idpdg = 99;
        id_epostopdg_(&_id[p], &idpdg);
        // Assignate entry variables (for precision concerns, issues if entry
        // variables used directly)
        double m = _m[p];
        double px = _px[p];
        double py = _py[p];
        double pz = _pz[p];
	
        // Calculate E
        double mt = sqrt(m * m + px * px + py * py);
        double E = mt * cosh(asinh(pz / mt));
        // Shift E and pz in the lab frame, if necessary
        if (*_dyframe != 0.0) {
          double E_lab = E * cosh(*_dyframe) - pz * sinh(*_dyframe);
          double pz_lab = pz * cosh(*_dyframe) - E * sinh(*_dyframe);
          E = E_lab;
          pz = pz_lab;
        }
        // Define the particle object and its mass
        HepMC3::GenParticlePtr particle = std::make_shared<HepMC3::GenParticle>(
            HepMC3::FourVector(px, py, pz, E), idpdg, 13);
        // Define the mass of the particle
        particle->set_generated_mass(m);

        // Signal if ID correspondance not found with 'id_epostopdg_'
        if (idpdg == 99) {
          std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°"
                    << evt.event_number() << " / vertex n°"
                    << particle->production_vertex()->id()
                    << ") :" << std::endl;
          std::cout << "\t Cannot find correspondent PDG_ID for decay child "
                       "particle with EPOS_ID = "
                    << _id[p] << " -> Recorded with PDG_ID = 99" << std::endl;
          std::cout << "     -> CHECK id_epostopdg_ (KW/ids.f) for this ID\n"
                    << std::endl;
        }

        // Mark it as 'read'
        used[p] = 1;
        // *** Attach the particle to the original vertex *** //
        v0->add_particle_out(particle);
      }
    }
    break;

  default:
    break;
  }

  // ***** Write the event ***** //
  writer->write_event(evt);
}

/**
 * @fn void closehepmc_()
 * @brief Function deleting the HepMC event object
 */
void closehepmc_() {
  std::cout << "***HepMC file close***" << std::endl;
  writer->close();
  delete writer;
}
