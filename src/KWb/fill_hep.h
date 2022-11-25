//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#ifdef WIN32
# ifdef CERNLIB_MSSTDCALL
#  define F77_UCASE
#  define type_of_call _stdcall
#  ifndef CERNLIB_QXCAPT
#    define CERNLIB_QXCAPT
#  endif
# else
#  define F77_LCASE
#  ifndef CERNLIB_QXNO_SC
#    define CERNLIB_QXNO_SC
#  endif
# endif
# define type_of_call  _stdcall
# define DEFCHARD   const char* , const int        
# define DEFCHARL          
# define PASSCHARD(string) string, strlen(string) 
# define PASSCHARL(string) 
#else
# define DEFCHARD     const char* 
# define DEFCHARL   , const int 
# define PASSCHARD(string) string 
# define PASSCHARL(string) , strlen(string) 
#endif
#ifdef CERNLIB_QXCAPT
#  define F77_NAME(name,NAME) NAME
#else
#  if defined(CERNLIB_QXNO_SC)
#    define F77_NAME(name,NAME) name
#  else
#    define F77_NAME(name,NAME) name##_
#  endif
#endif
#ifndef type_of_call
# define type_of_call
#endif


extern "C" void open_hepmc_(char* filename);
extern "C" void fillhepmc_(int *_iextree, int *_nevt, float *_eng, float *_dyframe, int *_iprojZ, int *_iprojA, int *_itargZ, int *_itargA, //*JJ 23/11
			int *_nhard, int *_ncoll, int *_npartproj, int *_nparttarg, int *_nspecp, int *_nspecn, int *_np, float *_bim, float *_sigtot, 
			int *_id, int *_ist, int *_ity, int *_ior, int *_jor, float *_zus, 
			float *_px, float *_py, float *_pz, float *_e, float *_x, float *_y, float *_z, float *_t); 
extern "C" void closehepmc_();

/// Call a Fortran subroutine
#define id_epostopdg F77_NAME(id_epostopdg,ID_EPOSTOPDG)
extern "C" {void type_of_call F77_NAME(id_epostopdg,ID_EPOSTOPDG)(int &idepos, int &idpdg);}   //**JJ
#define idtau F77_NAME(idtau,IDTAU)
extern "C" {void type_of_call F77_NAME(idtau,IDTAU)(int &id, float &p4, float &p5, float &tau);}   //**JJ

HepMC::IO_HEPEVT* hepevtio;
HepMC::IO_GenEvent* ascii_out;

