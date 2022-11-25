//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#pragma once

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <string.h>

using namespace std ;

// Q[] component indexes
#define T_ 0
#define X_ 1
#define Y_ 2
#define Z_ 3
#define NB_ 4
#define NQ_ 5
#define NS_ 6

#define PREDICT 0
#define CORRECT 1

#define C_1D 0
#define C_2D 1

const double C_PI=3.14159265358979312 ;


#define BAG_ASYMPT // necessary if you use HIRANO EoS !!!

const double gevtofm = 5.067728853 ;

#ifndef _DEBUG
//#define UI
#endif
