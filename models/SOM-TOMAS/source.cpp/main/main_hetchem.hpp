//=====================================================================================================
//                      THIS SUBROUTINE IMPLEMENTS THE SOM HETEROGENEOUS CHEMISTRY
//=====================================================================================================

void INTGR_HETCHEM(double *AECONC, double *params_in,double *OH,double *kOH_eff,double *Temp,double *Pres,double *dt, \
	           int *n_cells_in, int *CMAX_in, int *OMAX_in){
  
  // ACCOMMODATE INPUTS:
  // ==================================================================================================
  // NUMBER OF GRID CELLS:
  int n_cells = *n_cells_in;
  
  // GRID BOUNDS:
  int CMAX = *CMAX_in;
  int OMAX = *OMAX_in;
  
  // INTERNAL TIMESTEPS:
  //double ddt = 0.1;
  //int nsteps = (int)(*dt/ddt);

  // INTERNAL TIMESTEPS:
  int nsteps = 1; //(int)(*dt/ddt);
  double ddt = (*dt)/(double)nsteps;
  
  // THE SOM PARAMETERS:
  vector<double> params;
  params.push_back(params_in[0]);
  params.push_back(params_in[1]);
  params.push_back(params_in[2]);
  params.push_back(params_in[3]);
  params.push_back(params_in[4]);
  
  // INITIALIZE CELLS:
  // ==================================================================================================
  // INITIALIZE THE LISTS OF CELLS:
  vector<SOM_CELLS> list_cells;
  
  // INITIALIZE THE CELLS:
  for (int i = 1; i <= n_cells; i++){
    
    // CREATE CELL:
    SOM_CELLS cell(i,params,CMAX,OMAX,*Temp,*Pres,0);
    
    // READ kOH VALUE:
    cell.change_koh(*kOH_eff);
    
    // READ CONCENTRATION:
    cell.read_cgas(AECONC[i-1]);

    // ADD TO LIST:
    list_cells.push_back(cell);
  }
  
  // LOOP OVER TIMESTEPS AND SPECIES FOR CHEM.:  
  //=============================================
  for (int j = 1; j <= nsteps; j++){
    
    // CHEMICAL REACTION FOR CELLS:
    for (int i = 1; i <= n_cells; i++){
      
      SOM_CELLS &obj_cell = list_cells[i-1];

      obj_cell.step_chem(*OH,ddt,list_cells);
    }
    
    // SUM UP REACTIONS FOR CELLS:
    for (int i = 1; i <= n_cells; i++){
      
      SOM_CELLS &obj_cell = list_cells[i-1];

      obj_cell.sum_chem();
    }
  }
  
  // UPDATE THE CONC. ARRAY:
  //=============================================
  for (int i = 1; i <= n_cells; i++){
    
    SOM_CELLS &obj_cell = list_cells[i-1];

    AECONC[i-1] = obj_cell.Cgas;
  }
}

