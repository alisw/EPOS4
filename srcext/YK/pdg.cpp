//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//


#ifndef PARTICLE_PDG2
#include "pdg.h"
#endif

#include <iostream>
using std::cout;
using std::endl;

ParticlePDG2::ParticlePDG2() {
  fPDG   = kNonsensePDG;
  fMass  = -1.0;
  fWidth = 0.0;
  fNDecayChannels = 0;
};

ParticlePDG2::ParticlePDG2(Char_t *name, Int_t pdg, Double_t mass, Double_t width) {
  for(Int_t i=0; i<9; i++)
    if(*(name+i) != '\0') fName[i] = *(name+i);
    else break;
  fPDG   = pdg;
  fMass  = mass;
  fWidth = width;
  fNDecayChannels = 0;
};

ParticlePDG2::~ParticlePDG2() {
  for(Int_t i=0; i<fNDecayChannels; i++)
    delete fDecayChannels[i];
};

Double_t ParticlePDG2::GetFullBranching() {
  Double_t fullBranching = 0.0;
  for(Int_t i=0; i<fNDecayChannels; i++)
    fullBranching += fDecayChannels[i]->GetBranching();
  return fullBranching;
};

void ParticlePDG2::AddChannel(DecayChannel* channel) {
  if(channel->GetMotherPDG() != fPDG) {
    cout << "ERROR in ParticlePDG::AddChannel() : You try to add a channel which has a different mother PDG" << endl;
    return;
  }
  fDecayChannels[fNDecayChannels] = new DecayChannel() ;
  fDecayChannels[fNDecayChannels]->SetMotherPDG(channel->GetMotherPDG());
  fDecayChannels[fNDecayChannels]->SetBranching(channel->GetBranching());
  fDecayChannels[fNDecayChannels]->SetDaughters(channel->GetDaughters(), channel->GetNDaughters());
  fNDecayChannels++;
};
