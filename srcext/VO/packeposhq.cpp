//
// This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
// (See COPYING file for the text of the licence)
//

extern "C" {
  void hqini_();
}


/**
 * Initialize heavy quark part
 *
 */
void initializeHeavyQuarkPart(){  
  hqini_();
}
