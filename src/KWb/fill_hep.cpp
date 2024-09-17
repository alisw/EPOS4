//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/HepMCDefs.h"
#include "HepMC/IO_HEPEVT.h"//*JJ
#include "HepMC/HEPEVT_Wrapper.h" //*JJ
#include "HepMC/IO_AsciiParticles.h"//*JJ
#include <iostream>
#include <cstdlib>
#include <math.h> //*JJ
#include <stdlib.h> //**JJ
#include <stdio.h> //**JJ

#include "fill_hep.h"

HepMC::IO_GenEvent* writer;
bool is_std = false;

/// ***** Function creating the HepMC event object ***** ///
void open_hepmc_(char* filename){
    if (filename == "/dev/stdout")
        is_std = true;
	writer = new HepMC::IO_GenEvent(filename);
    if(!writer){std::cout << " FAILED!!!! "<< std::endl; }
    std::cout << "*** HepMC 2.06.09 ***" << std::endl;  
    std::cout << "*** " << filename << " open***" << std::endl; 
}

/// *** Function assigning the mass to the colliding ions *** ///
double get_mass(int A, int Z){
	double m = 0;
//U   : 
	if(A==238 && Z==92) m=221.742894;
//Pb : 
	if(A==208 && Z==82) m=193.729016;
//Au : 
	if(A==197 && Z==79) m=183.473190;
//Xe : 
	if(A==129 && Z==54) m=120.074038;
//Cu : 
	if(A==63 && Z==29)  m=58.6185461;
//d   : 
	if(A==2 && Z==1)    m=1.87612384;
//p  :
	if(A==1 && Z==1)    m=0.93827208;
return m;	
}

/// ***** Function filling the HepMC event ***** ///
void fillhepmc_(int *_iextree, int *_nevt, float *_eng, float *_dyframe, int *_iprojZ, int *_iprojA, int *_itargZ, int *_itargA,
		int *_nhard, int *_ncoll, int *_npartproj, int *_nparttarg, int *_nspecp, int *_nspecn, int *_np, float *_bim, float *_sigtot, 
		int *_id, int *_ist, int *_ity, int *_ior, int *_jor, float *_zus, float *_px, float *_py, float *_pz,  float *_e, 
		float *_x, float *_y, float *_z, float *_t){

    /// *********************************** ///     
	 /// ****** EVENT INITIALISATION ******* ///
  /// *********************************** ///
  
  /// *** Event definition *** ///
	HepMC::GenEvent evt; //create event
	evt.use_units(HepMC::Units::GEV, HepMC::Units::MM); //units
	evt.set_event_number(*_nevt); //event number
	
	/// *** Cross-section *** ///
	HepMC::GenCrossSection cs; //generate the cross-section object
	if (*_iextree == 1){ //provide value only if new tree used (from EPOS 4)
		cs.set_cross_section(*_sigtot * 1e9); //mb (EPOS) -> pb (HepMC)
	} else {
		cs.set_cross_section(0.0);
	}
	evt.set_cross_section(cs);
	
	/// *** Beam information *** ///
	int id_projectile, id_target;
	// A,Z
	int A_proj = *_iprojA;
	int Z_proj = *_iprojZ;
	int A_targ = *_itargA;
	int Z_targ = *_itargZ;
	// ID
	if(Z_proj == 1 && A_proj== 1) {
		id_projectile = 2212;
	} else { 
		id_projectile = 1000000000 + Z_proj*10000 + A_proj*10;
	}
	if(Z_targ == 1 && A_targ == 1) {
		id_target = 2212;
	} else {
		id_target = 1000000000 + Z_targ*10000 + A_targ*10;
	}
		
	/// *** Heavy ion information *** ///
	if(A_proj > 1 || A_targ > 1){
		HepMC::HeavyIon hi;
		hi.set_impact_parameter(*_bim); //impact parameter
		// Provide following info only if new tree used (from EPOS 4)
		if (*_iextree == 1){ 
			hi.set_Ncoll_hard(*_nhard); //number of hard interactions
			hi.set_Ncoll(*_ncoll); //number of Glauber nucleon collisions
			hi.set_Npart_proj(*_npartproj); //number of participants in projectile
			hi.set_Npart_targ(*_nparttarg); //number of participants in target
			hi.set_spectator_protons(*_nspecp); //number of protons spectators
			hi.set_spectator_neutrons(*_nspecn); //number of neutrons spectators
		}
		evt.set_heavy_ion(hi);
	}

	/// *** Defining the beam particles *** ///
	// Mass
	double mass_beam1 = get_mass(A_targ ,Z_targ);
	double mass_beam2 = get_mass(A_proj ,Z_proj);
	
	// Energy
	double energy_beam1 = *_eng * 0.5 * A_targ ;
	double energy_beam2 = *_eng * 0.5 * A_proj ;
	
	if(*_dyframe != 0.0){
		if(*_nevt == 1) std::cout << "  Particles shifted from CoM to lab with dy=" << *_dyframe << std::endl;  
		energy_beam1 = cosh(*_dyframe) * energy_beam1 - sinh(*_dyframe) * (-1)*sqrt(energy_beam1 * energy_beam1 - mass_beam1 * mass_beam1);
		energy_beam2 = cosh(*_dyframe) * energy_beam2 - sinh(*_dyframe) * sqrt(energy_beam2 * energy_beam2 - mass_beam2 * mass_beam2);
	}
	
	// Create the corresponding particle objects
	HepMC::GenParticle* gp_beam1 = new HepMC::GenParticle(HepMC::FourVector(0,0,(-1)*sqrt(energy_beam1 * energy_beam1 - mass_beam1 * mass_beam1), energy_beam1), id_target, 4);
	HepMC::GenParticle* gp_beam2 = new HepMC::GenParticle(HepMC::FourVector(0,0,sqrt(energy_beam2 * energy_beam2 - mass_beam2 * mass_beam2), energy_beam2), id_projectile, 4);
	
	// Assign the masses
	gp_beam1->set_generated_mass(mass_beam1);
	gp_beam2->set_generated_mass(mass_beam2);
	
	// Attach them to the event object
	evt.set_beam_particles(gp_beam1, gp_beam2);

	/// ***** Zero Vertex ***** ///
	HepMC::GenVertex* v0 = new HepMC::GenVertex(HepMC::FourVector(0, 0, 0, 0) ); 
	evt.add_vertex(v0); 
	
	/// *** Adding colliding nuclei to vertex 0 *** ///
	v0 -> add_particle_in(gp_beam1);
	v0 -> add_particle_in(gp_beam2);
	
	/// *** Set barcodes for beams *** ///
	evt.set_signal_process_id(v0->id());
	evt.set_signal_process_vertex(v0);
	
	  /// ******************************** ///
	 /// *** DEFINE GENERAL VARIABLES *** ///
	/// ******************************** ///
	// Total number of particles
	int np = *_np;
	
	// PDG ID variable
	int idpdg;
	
	// Particle lifetime
	float tau;
	
	// Define the "used" list initialised with -1 :
	// - "-1" = never used
	// -  "0" = used as vtx output but NOT final-state
	// -  "1" = final-state OR used as vtx output + input
	int used[np];
	for (int i = 0; i < np; i++) used[i] = -1;
	
	// Define a vector of pointers which will contain decay children that will decay again
	// (so avoid declaring twice the same particle)
	std::vector<HepMC::GenParticle*> temp_ptl;
	
	// Define which particle selection is applied (1 : only FS particles || 2 : FS particles + parents & vertex for weak & EM decayed particles)
	int which_table = 2 ;
	
	// Limit lifetime to distinguish strong from weak/EM decay
	float tau_wkEM = 1.e-19;
	
	// List of particles ID with tau < tau_wkEM but for which decay record is needed 
	// Contains : [J/Psi, Psi(2S), Psi(3770), Upsilon, Upsilon(2S), Upsilon(3S), Upsilon(4S)]
	const int num_ex = 7;
	int excptn_list [num_ex] = {441, 447, 448, 551, 556, 557, 558};
	
	// Default energy and mass to have gamma = 1 in idtau(..)
	float p4 = 1.0, p5 = 1.0;
	
	/// ***** Selection of the EPOS particles written in the HepMC event ***** ///
	switch(which_table){
		
		  /// ************************************************ ///
		 /// ***** Final-state particles only (ist = 0) ***** ///
		/// ************************************************ ///
		case 1:
			
			/// *** Loop over the EPOS particle list *** ///
			for(int p = 0; p < np; p++){
				if (_ist[p] == 0 && used[p] == -1){
					// Translate ID
					idpdg = 99;
					id_epostopdg_(_id[p],idpdg);
					// Assignate entry variables (for precision concerns, issues if entry variables used directly)
					double m = _e[p];
					double px = _px[p];
					double py = _py[p];
					double pz = _pz[p];
					// Calculate E
					double mt = sqrt(m*m + px*px + py*py);
					double E = mt * cosh( asinh(pz/mt) );
					// Shift E and pz in the lab frame, if necessary
					if(*_dyframe != 0.0){
						double E_lab  = E  * cosh(*_dyframe) - pz * sinh(*_dyframe);
						double pz_lab = pz * cosh(*_dyframe) - E  * sinh(*_dyframe);
						E  = E_lab ;
						pz = pz_lab ;
					}
					// Define the particle object and its mass
					HepMC::GenParticle* particle = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), idpdg, 1);
					// Define the mass of the particle
					particle -> set_generated_mass(m);
					
					// Signal if ID correspondance not found with 'id_epostopdg'
					if (idpdg == 99) {
						std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°" << evt.event_number() << " / vertex n°" << particle -> production_vertex() -> barcode() << ") :" << std::endl ;
						std::cout << "\t Cannot find correspondent PDG_ID for decay child particle with EPOS_ID = " << _id[p] << " -> Recorded with PDG_ID = 99" << std::endl ;
						std::cout << "     -> CHECK id_epostopdg (KW/ids.f) for this ID\n" << std::endl ;
					}

					// Mark it as 'read'
					used[p] = 1;
					/// *** Attach the particle to the original vertex *** ///
					v0 -> add_particle_out(particle);
				}
			}
			
			break;
		
		  /// ****************************************************************** ///
		 /// ***** Final-state particles + decays (weak, EM & exceptions) ***** ///
		/// ****************************************************************** ///
		case 2:
			
			/// *** Loop over the EPOS particle list *** ///
			for(int p = 0; p < np; p++){

				/// ***** For weak or EM decayed particles ***** /// 
				if (_ist[p] == 1 && used[p] != 1) {
					
					/// *** Determines the lifetime*** ///
					idtau(_id[p],p4,p5,tau);
					//Calculate tau in seconds from c*tau in fm
					tau = tau * 1.e-15 / 3.0e8;
					
					/// *** Determines the decay type from parent's lifetime *** ///
					// To mark if decay is strong or weak/EM
					bool weak_EM = false;
					// Evaluate the decay type (EXCEPT FOR INTERACTING NUCLEONS)
					if (tau > tau_wkEM && p > A_targ+A_proj-1) weak_EM = true;

					/// *** Record the decay if particle in exceptions list *** ///
					for(int e = 0; e < num_ex; e++){
						if(abs(_id[p]) == excptn_list[e]){
							weak_EM = true;
						}
					}
					
					/// *** STRONG DECAY + INTERACTING NUCLEON + used[p] != 0 -> ignored *** ///
					if(weak_EM == false && used[p] == -1){
						
						// Mark parent as 'read'
						used[p] = 1;
						continue;

					/// *** WEAK / EM DECAY + EXCEPTIONS -> recorded *** ///
					} else {
						/// *** Parent's information *** ///
						// Translate ID
						idpdg = 99;
						id_epostopdg_(_id[p],idpdg);
						// Assignate entry variables (for precision concerns, issues if entry variables used directly)
						double m = _e[p];
						double px = _px[p];
						double py = _py[p];
						double pz = _pz[p];
						// Calculate E
						double mt = sqrt(m*m + px*px + py*py);
						double E = mt * cosh( asinh(pz/mt) );
						// Shift E and pz in the lab frame, if necessary
						if(*_dyframe != 0.0){
							double E_lab  = E  * cosh(*_dyframe) - pz * sinh(*_dyframe);
							double pz_lab = pz * cosh(*_dyframe) - E  * sinh(*_dyframe);
							E  = E_lab ;
							pz = pz_lab ;
						}
						// Define the particle object and its mass
						HepMC::GenParticle* parent = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), idpdg, 2);
						
						// Define the mass of the particle
						parent -> set_generated_mass(m);
						
						// If parent is NOT originating from a decay
						if (used[p] == -1) {
							
							/// *** Attach the parent as a child of the original vertex if not from a decay *** ///
							v0 -> add_particle_out(parent);
							
						// If parent is originating from a decay
						} else {

							// Search for it in the temporary storing list 'temp_ptl'
							for(std::size_t t = 0; t < temp_ptl.size(); t++){
								// Test if parent is the particle checked in 'temp_ptl' (compare mass, energy and px)
								if(m == temp_ptl[t]->generated_mass() && E == temp_ptl[t]->momentum().e() && px == temp_ptl[t]->momentum().px()){
									// Test if ID is the same
									if (idpdg != temp_ptl[t]->pdg_id()) {
										std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°" << evt.event_number() << ") :" << std::endl ;
										std::cout << "     Found correspondence between current particle and a decay child (PDG_ID = " << temp_ptl[t]->pdg_id() <<"), but ID doesn't match" << std::endl ;
										std::cout << "\t Current particle : EPOS_ID = " << _id[p] << "  /  corresponding PDG_ID = " << idpdg << " (from id_epostopdg)" << std::endl ;
										std::cout << "     -> Recorded as the same particle anyway, but CHECK id_epostopdg (KW/ids.f) for this ID\n" << std::endl ;
									}
									// Assign 'parent' to its already created corresponding object
									parent = temp_ptl[t];
									// Remove the object from the 'temp_ptl' list
									temp_ptl.erase(temp_ptl.begin()+t);
								}
							}
							
						}
						
						/// *** Search for the children *** ///
						// Marker used to indicate if vertex already defined
						bool vertex = false;
						// Marker used to indicate when to stop reading
						bool family = false;
						
						// Create the pointer for secondary vertex
						HepMC::GenVertex* v = new HepMC::GenVertex();
						
						// Loop over the other particles
						for(int c = p+1; c < np; c++){
							
							// Check if child from particle p
							if(_ior[c] == p+1){
								
								// Indicate that a child is found
								family = true;
								
								/// *** Secondary vertex *** ///
								if(vertex == false){
									// Convert units from fm (EPOS) to mm (HepMC)
									double x = _x[c] * 1.0e-12 ;
									double y = _y[c] * 1.0e-12 ;
									double z = _z[c] * 1.0e-12 ;
									double t = _t[c] ;
									// Define the vertex
									v -> set_position(HepMC::FourVector(x, y, z, t));
									// Attach it to the event
									evt.add_vertex(v);
									
									// Mark vertex as defined
									vertex = true;
									
									/// *** Attach the parent to the secondary vertex *** ///
									v -> add_particle_in(parent);
								}

								/// *** Child's information *** ///
								// Translate ID
								idpdg = 99;
								id_epostopdg_(_id[c],idpdg);
								// Assign entry variables (for precision concerns, issues if entry variables used directly)
								double m = _e[c];
								double px = _px[c];
								double py = _py[c];
								double pz = _pz[c];
								// Calculate E
								double mt = sqrt(m*m + px*px + py*py);
								double E = mt * cosh( asinh(pz/mt) );
								// Shift E and pz in the lab frame, if necessary
								if(*_dyframe != 0.0){
									double E_lab  = E  * cosh(*_dyframe) - pz * sinh(*_dyframe);
									double pz_lab = pz * cosh(*_dyframe) - E  * sinh(*_dyframe);
									E  = E_lab ;
									pz = pz_lab ;
								}
								// Define the particle object for child
								HepMC::GenParticle* child;
								
								// K0(b) (EXCEPTION) : do not record it but directly its K0s/L child
								int c_K0 = c; //record the index in order to come back to it at the end
								bool K0 = false; //to indicate if child is a K0 that has been ignored
								
								if(abs(idpdg) == 311) {
									// Indicate that a K0(b) has been found
									K0 = true;
									// Search for the child K0s/L 
									for(int K0sl = c+1; K0sl < np; K0sl++) {
										if(_ior[K0sl] == c+1){
											// Mark the K0(b) as 'used'
											used[c] = 1;
											// Change index to the K0s/L one (for particle information)
											c = K0sl; 
											// Translate ID of the K0s/L
											id_epostopdg_(_id[K0sl],idpdg);
											// Break the loop
											break;
										}
									}
								}

								if (_ist[c] == 0) {
									// If FS particle : set hep_status = 1 
									child = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), idpdg, 1);
									// Mark child as 'final-state'
									used[c] = 1;
								} else {
									// If not FS particle (to decay later) : keep the object until used as parent + set hep_status = 2 
									temp_ptl.push_back(new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), idpdg, 2));
									// Assign 'child' to its corresponding stored object to be re-used
									child = temp_ptl[temp_ptl.size()-1];
									// Mark child as 'used as input but NOT final-state'
									used[c] = 0;
								}
								// Define the mass of the child
								child -> set_generated_mass(_e[c]);
								
								/// *** Attach the child to the secondary vertex *** ///
								v -> add_particle_out(child);
								
								// K0(b) (EXCEPTION)
								if(K0 == true) {
									// Re-assign 'c' counter to the originally found K0(b) index
									c = c_K0;
								}
								
								// Signal if ID correspondance not found with 'id_epostopdg'
								if (idpdg == 99) {
									std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°" << evt.event_number() << " / vertex n°" << child -> production_vertex() -> barcode() << ") :" << std::endl ;
									std::cout << "\t Cannot find correspondent PDG_ID for decay child particle with EPOS_ID = " << _id[c] << " -> Recorded with PDG_ID = 99" << std::endl ;
									std::cout << "     -> CHECK id_epostopdg (KW/ids.f) for this ID\n" << std::endl ;
								}

							// Check if children have been found already
							} else if(_ior[c] != p+1 && family == true) {
									// Break the loop if children already found 
									// (uses the fact that children from same parent always follow each other in the list)
									break;
								
							// Read the next particle if not a child
							} else {
								continue;
							}
						}
						
						/// *** No children found : ERROR MESSAGE *** ///
						if(family == false) {
							std::cout << "\n /!\\ NO CHILDREN FOUND DESPITE ist=1 FOR PARTICLE n°" << p << " (ID = " << _id[p] << ") \n" << std::endl;
							// Set hep_status = 1 then (FS particle)
							parent -> set_status(12);
						}
						
						// Mark parent as 'FS / used for output + input'
						used[p] = 1;
					}
					
				/// ***** For other final-state particles (not from weak or EM decay) ***** ///
				} else if (_ist[p] == 0 && used[p] == -1){
					// Translate ID
					idpdg = 99;
					id_epostopdg_(_id[p],idpdg);
					// Assign entry variables (for precision concerns, issues if entry variables used directly)
					double m = _e[p];
					double px = _px[p];
					double py = _py[p];
					double pz = _pz[p];
					// Calculate E
					double mt = sqrt(m*m + px*px + py*py);
					double E = mt * cosh( asinh(pz/mt) );
					// Shift E and pz in the lab frame, if necessary
					if(*_dyframe != 0.0){
						double E_lab  = E  * cosh(*_dyframe) - pz * sinh(*_dyframe);
						double pz_lab = pz * cosh(*_dyframe) - E  * sinh(*_dyframe);
						E  = E_lab ;
						pz = pz_lab ;
					}
					// Define the particle object and its mass
					HepMC::GenParticle* child = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), idpdg, 1);
					// Define the mass of the particle
					child -> set_generated_mass(m);
					// Mark it as 'read'
					used[p] = true;
					/// *** Attach the particle to the original vertex *** ///
					v0 -> add_particle_out(child);
					
					// Signal if ID correspondance not found with 'id_epostopdg'
					if (idpdg == 99) {
						std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°" << evt.event_number() << " / vertex n°" << child -> production_vertex() -> barcode() << ") :" << std::endl ;
						std::cout << "\t Cannot find correspondent PDG_ID for decay child particle with EPOS_ID = " << _id[p] << " -> Recorded with PDG_ID = 99" << std::endl ;
						std::cout << "     -> CHECK id_epostopdg (KW/ids.f) for this ID\n" << std::endl ;
					}
				}
			}
			
			break;
					
		  /// ******************************************** ///
		 /// ***** Hadrons before hadronic cascades ***** ///
		/// ******************************************** ///
		case 3:
			
			/// *** Loop over the EPOS particle list *** ///
			for(int p = 0; p < np; p++){
				if (_ist[p] == 3 && used[p] == -1){
					// Translate ID
					idpdg = 99;
					id_epostopdg_(_id[p],idpdg);
					// Assignate entry variables (for precision concerns, issues if entry variables used directly)
					double m = _e[p];
					double px = _px[p];
					double py = _py[p];
					double pz = _pz[p];
					// Calculate E
					double mt = sqrt(m*m + px*px + py*py);
					double E = mt * cosh( asinh(pz/mt) );
					// Shift E and pz in the lab frame, if necessary
					if(*_dyframe != 0.0){
						double E_lab  = E  * cosh(*_dyframe) - pz * sinh(*_dyframe);
						double pz_lab = pz * cosh(*_dyframe) - E  * sinh(*_dyframe);
						E  = E_lab ;
						pz = pz_lab ;
					}
					// Define the particle object and its mass
					HepMC::GenParticle* particle = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), idpdg, 1);
					// Define the mass of the particle
					particle -> set_generated_mass(m);
					
					// Signal if ID correspondance not found with 'id_epostopdg'
					if (idpdg == 99) {
						std::cout << "\n /!\\ WARNING (from 'fill_hep.cpp', event n°" << evt.event_number() << " / vertex n°" << particle -> production_vertex() -> barcode() << ") :" << std::endl ;
						std::cout << "\t Cannot find correspondent PDG_ID for decay child particle with EPOS_ID = " << _id[p] << " -> Recorded with PDG_ID = 99" << std::endl ;
						std::cout << "     -> CHECK id_epostopdg (KW/ids.f) for this ID\n" << std::endl ;
					}

					// Mark it as 'read'
					used[p] = 1;
					/// *** Attach the particle to the original vertex *** ///
					v0 -> add_particle_out(particle);
				}
			}
			
			break;
					
		default:
			break;
	}

  /// ***** Write the event ***** ///
	writer -> write_event(&evt);
  
}

/// ***** Function deleting the HepMC event object ***** ///
void closehepmc_(){
	if (is_std == false)
	 std::cout << "***hepmc close***" << std::endl;  
   delete writer;
}

