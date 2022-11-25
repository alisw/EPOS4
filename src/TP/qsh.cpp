//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

extern "C" {
  void gammaa_();
}


/**
 * Initialize electron proton part
 * fake e-A collision with gamma(pi0)-A
 *
 */
void initializeElectronProtonPart(){
  gammaa_();
}                              
