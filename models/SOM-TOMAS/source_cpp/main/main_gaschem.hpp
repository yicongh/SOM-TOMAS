//=====================================================================================================
//                       THIS SUBROUTINE IMPLEMENTS THE SOM GAS-PHASE CHEMISTRY
//=====================================================================================================

void INTGR_GASCHEM(double *SAPRCGC, double *params_in,double *OH, double *Temp, double *Pres, double *dt, \
	           int *n_precs_in, int *n_cells_in, int *i1_precs_in, int *i1_cells_in, \
	           int *CNO_prec_in, int *ONO_prec_in, double *kOH_prec_in, \
	           int *CMAX_in, int *OMAX_in, int *G1ONLY){
  
  // ACCOMMODATE INPUTS:
  // ==================================================================================================
  int n_precs = *n_precs_in;
  int n_cells = *n_cells_in;
  
  int i1_precs = *i1_precs_in;
  int i1_cells = *i1_cells_in;

  vector<int>    CNO_prec;
  vector<int>    ONO_prec;
  vector<double> kOH_prec;
  
  for (int i = 1; i <= n_precs; i++){
    CNO_prec.push_back(CNO_prec_in[i-1]);
    ONO_prec.push_back(ONO_prec_in[i-1]);
    kOH_prec.push_back(kOH_prec_in[i-1]);
  }

  int CMAX = *CMAX_in;
  int OMAX = *OMAX_in;

  // INTERNAL TIMESTEPS:
  int nsteps = 10; //(int)(*dt/ddt);
  double ddt = (*dt)/(double)nsteps;

  // THE SOM PARAMETERS:
  vector<double> params;
  params.push_back(params_in[0]);
  params.push_back(params_in[1]);
  params.push_back(params_in[2]);
  params.push_back(params_in[3]);
  params.push_back(params_in[4]);

  // INITIALIZE CELLS AND PRECURSORS:
  // ==================================================================================================
  // INITIALIZE THE LISTS OF CELLS AND PRECURSORS:
  vector<SOM_PRECS> list_precs;
  vector<SOM_CELLS> list_cells;

  // INITIALIZE THE PRECURSORS:
  for (int i = 1; i <= n_precs; i++){
   
    int CNO = CNO_prec[i-1];
    int ONO = ONO_prec[i-1];
    
    // CREATE PRECURSOR:
    SOM_PRECS prec(CNO,ONO,params,CMAX,OMAX,*Temp,*Pres);
    
    // ADD TO LIST:
    list_precs.push_back(prec);
  }
  
  // INITIALIZE THE CELLS:
  for (int i = 1; i <= n_cells; i++){
    
    // CREATE CELL:
    SOM_CELLS cell(i,params,CMAX,OMAX,*Temp,*Pres,*G1ONLY);
    
    // ADD TO LIST:
    list_cells.push_back(cell);
  }

  //exit(0);
  
  // READ kOH AND GAS CONCENTRATION FOR PRECURSORS:
  for (int i = 1; i <= n_precs; i++){
    
    // SELECT PRECURSOR:
    SOM_PRECS &obj_prec = list_precs[i-1];

    // READ kOH:
    obj_prec.read_koh(kOH_prec[i-1]);

    // READ CONCENTRATION:
    obj_prec.read_cgas(SAPRCGC[i1_precs+i-1]);
    
  }
  
  // READ CONCENTRATION FOR CELLS:
  for (int i = 1; i <= n_cells; i++){

    // SELECT CELL:
    SOM_CELLS &obj_cell = list_cells[i-1];

    // READ CONCENTRAION:
    obj_cell.read_cgas(SAPRCGC[i1_cells+i-1]);
    
  }
  
  // STEP FOR GAS-PHASE CHEMISTRY:
  // ==================================================================================================
  for (int j = 1; j <= nsteps; j++){
  
    // CHEMICAL REACTION FOR PRECURSORS:
    for (int i = 1; i <= n_precs; i++){
      
      SOM_PRECS &obj_prec = list_precs[i-1];
    
      obj_prec.step_chem(*OH,ddt,list_cells);
  
    }
    
    // CHEMICAL REACTION FOR CELLS:
    for (int i = 1; i <= n_cells; i++){
      
      SOM_CELLS &obj_cell = list_cells[i-1];

      obj_cell.step_chem(*OH,ddt,list_cells);
   
    }
    
    // SUM UP REACTIONS FOR PRECURSORS:
    for (int i = 1; i <= n_precs; i++){
      
      SOM_PRECS &obj_prec = list_precs[i-1];

      obj_prec.sum_chem();
    
    }
    
    // SUM UP REACTIONS FOR CELLS:
    for (int i = 1; i <= n_cells; i++){

      SOM_CELLS &obj_cell = list_cells[i-1];

      obj_cell.sum_chem();

    }
  }
  
  // UPDATE SAPRCGC BASED ON RESULTS:
  // ==================================================================================================
  for (int i = 1; i <= n_precs; i++){
    
    SOM_PRECS &obj_prec = list_precs[i-1];

    SAPRCGC[i1_precs+i-1] = obj_prec.Cgas;
 
  }

  for (int i = 1; i <= n_cells; i++){
    
    SOM_CELLS &obj_cell = list_cells[i-1];

    SAPRCGC[i1_cells+i-1] = obj_cell.Cgas;
  }
    
}


