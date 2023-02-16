// =========================================================================================================
//                THIS IS THE INTERFACE FOR THE C++ MODULE FOR GAS-PHASE AND HET. CHEMISTRY
// =========================================================================================================

#include "main_header.h"
#include "main_gaschem.hpp"
#include "main_hetchem.hpp"

// GAS-PHASE CHEMISTRY:
// =========================================================================================================
extern"C" double gas_chem_(double *saprcgc,double *params,double *OH,double *Temp,double *Pres,double *dt, \
			   int *n_precs_in, int *n_cells_in, int *i1_precs_in, int *i1_cells_in,           \
			   int *CNO_prec_in, int *ONO_prec_in, double *kOH_prec_in,                        \
			   int *CMAX_in, int *OMAX_in, int *G1ONLY){
  
  INTGR_GASCHEM(saprcgc,params,OH,Temp,Pres,dt,		  \
	        n_precs_in,n_cells_in,i1_precs_in,i1_cells_in, \
	        CNO_prec_in,ONO_prec_in,kOH_prec_in,           \
	        CMAX_in,OMAX_in,G1ONLY);
  
  return 0.;
}

// HETEROGENEOUS CHEMISTRY:
// =========================================================================================================
extern"C" double het_chem_(double *aeconc,double *params,double *OH,double *kOH_eff, \
			   double *Temp, double *Pres,double *dt,	             \
			   int *n_cells_in, int *CMAX_in, int *OMAX_in){

  INTGR_HETCHEM(aeconc,params,OH,kOH_eff,Temp,Pres,dt, \
	        n_cells_in,CMAX_in,OMAX_in);
  
  return 0.;
}
