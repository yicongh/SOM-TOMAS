//===============================================
// This file defines the class for the precursors
//===============================================

#pragma once

// Define the class:
class SOM_PRECS{

public:
  // Precursor carbon and oxygen number:
  int CNO,ONO; 
  int CMAX,OMAX;
  double T,P;
  vector<double> params;
  
  SOM_PRECS(int CNO_in,int ONO_in,vector<double> params_in, \
	    int CMAX_in,int OMAX_in,                        \
	    double Temp_in, double Pres_in){
    
    CNO    = CNO_in;
    ONO    = ONO_in;
    CMAX   = CMAX_in;
    OMAX   = OMAX_in;
    T      = Temp_in;
    P      = Pres_in;
    params = params_in;

    auto_0();
    auto_1();
    auto_2();
    auto_3();
    auto_4();
    auto_5();
    auto_6();
    
  };

  // Extract environment variables:
  //int CMAX,OMAX;

  void auto_0(){

    //CMAX = env::CMAX;
    //OMAX = env::OMAX;

    //P = env::P;
    //T = env::T;
  }

  // Find ID of precursor in grid (if ONO > 0):
  int ID;
  void auto_1(){
    
    if (ONO == 0){
      ID = calc_id(CNO,ONO+1,CMAX,OMAX);
    } else {
      ID = calc_id(CNO,ONO,CMAX,OMAX);
    }
    
  };
  
  // Molecular weight:
  double MW;
  void auto_2(){
    MW = double(CNO)*12.0107 + double(ONO)*15.999 + ((double(CNO)*2) + 2 - double(ONO))*1.00794;
  };
  
  // Find the number of dependants of functionalization & fragmentation:
  int n_func,n_frag,x_frag;
  void auto_3(){
    if (ONO > 0){
      n_func             = find_nfunc(ID,CNO,ONO,CMAX,OMAX);
      //tie(n_frag,x_frag) = find_nfrag(ID,CNO,ONO,CMAX,OMAX);
      n_frag  = find_nfrag(ID,CNO,ONO,CMAX,OMAX);
    } else {
      n_func = min(4,CNO*2);
      n_frag = 0;
    }
  };
  
  // For functionalization:
  vector<int> depend_func;
  void auto_4(){
    if (ONO > 0){
      depend_func = find_depend_func(ID,n_func);
    } else {
      for (int i = 0;i <= n_func - 1;i++){
	depend_func.push_back(ID + i);
      }
    }
  };

  // For fragmentation:
  vector<int> depend_frag;
  void auto_5(){
    if (ONO > 0){
      depend_frag = find_depend_frag(ID,CNO,ONO,CMAX,OMAX,n_frag);
    } else {
      ;
    }
  };

  // Extract the parameters:
  double pfrag,pfunc;
  vector<double> pfs;
  void auto_6(){

    double mfrag = params[0];
    pfrag = pow(double(ONO)/double(CNO),mfrag);
    pfunc = 1. - pfrag;

    pfs.push_back(params[1]);
    pfs.push_back(params[2]);
    pfs.push_back(params[3]);
    pfs.push_back(params[4]);
  }

  // Read reactivity with OH:
  double kOH,kOH_ppm;
  void read_koh(double kOH_in){

    kOH = kOH_in;

    double Na = 6.02e23; // [molecules mol-1]
    double R  = 8.314;   // [J mol-1 K-1]

    kOH_ppm = kOH*(7.3395e15)/300.; //*(1e-6)*Na*(101325./(R*300.))*(1e-6); // [ppm-1 s-1]
  };

  // Assign kOH to corresponding cell:
  void assign_koh(vector<SOM_CELLS> &list_cells){

    if (ONO > 0){
      int ID_assign = calc_id(CNO,ONO,CMAX,OMAX);

      // Find the cell to assign to:
      SOM_CELLS &cell_assign = list_cells[ID_assign-1];

      // Now assign kOH to cell:
      cell_assign.change_koh(kOH_ppm);
    } else {
      ;
    }
  }
  
  
  // Read gas-phase concentration at initial time:
  double Cgas;
  void read_cgas(double Cgas_in){

    Cgas = Cgas_in;
  }
  
  // Initiate add and del variables:
  double to_add = 0.;
  double to_del = 0.;

  // Gas-phase reaction:
  void step_chem(double OH_in,double dt,vector<SOM_CELLS> &list_cells){
    // Convert OH concentration to ppm:
    double Na = 6.02e23; // [molecules mol-1]
    double R  = 8.314;   // [J mol-1 K-1]

    double OH;
    OH = OH_in;
    
    // The amount of species reacted:
    double k1,k2,k3,k4;

    k1 = kOH_ppm*OH*Cgas*(dt);
    k2 = kOH_ppm*OH*(Cgas+0.5*k1)*(dt);
    k3 = kOH_ppm*OH*(Cgas+0.5*k2)*(dt);
    k4 = kOH_ppm*OH*(Cgas+k3)*(dt);
    
    to_del = to_del + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
    
    // Add to dependants of functionalization:
    int ID_add;
    
    for (int i = 0;i <= n_func-1;i++){

      ID_add = depend_func[i];

      // Find the cell to add to:
      SOM_CELLS &cell_add = list_cells[ID_add-1];

      // Now add to cell:
      cell_add.to_add = cell_add.to_add + to_del*pfunc*pfs[i];
    }

    // Add to dependants of fragmentation:
    for (int i = 0;i <= n_frag-1;i++){

      ID_add = depend_frag[i];

      // Find the cell to add to:
      SOM_CELLS &cell_add = list_cells[ID_add-1];

      // Now add to cell:
      double fac = 1./double(n_frag);
      fac = roundf(fac*1000)/1000;
      
      cell_add.to_add = cell_add.to_add + to_del*pfrag*fac;
    }
 
  };

  // Sum up addition and reduction:
  void sum_chem(){

    Cgas = Cgas + to_add - to_del;
    
    to_add = 0.;
    to_del = 0.;
  }

};

