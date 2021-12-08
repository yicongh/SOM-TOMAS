//===================================================
// THIS FILE DEFINES THE CLASS FOR THE SOM GRID CELLS
//===================================================

#pragma once

// DEFINE THE CLASS:
//===================================================
class SOM_CELLS{

public:
  
  // INSTANCE INITILIZATION:
  //=================================================
  // SOM GRID ID AND,
  // SOM PARAMETERS:
  int ID;
  vector<double> params;
  
  // MAX. CARBON AND OXYGEN NUMBERS AND
  int CMAX,OMAX;

  int G1ONLY;
  
  // CHAMBER TEMPERATURE [K] AND
  // PRESSURE [Pa]:
  double T,P;

  SOM_CELLS(int ID_in,vector<double> params_in, \
	    int CMAX_in,int OMAX_in,            \
	    double Temp_in,double Pres_in,int G1ONLY_IN){
   
    ID     = ID_in;
    CMAX   = CMAX_in;
    OMAX   = OMAX_in;
    T      = Temp_in;
    P      = Pres_in;
    params = params_in;
    G1ONLY = G1ONLY_IN;

    // AUTOMATIC ROUTINES:
    auto_0();
    auto_1();
    auto_2();
    auto_3();
    auto_4();
    auto_5();
    auto_6();
    auto_7();
    
  };

  // EXTRACT ENVIRONMENT VARIABLES:
  //=================================================
  // PRESSURE AND
  // TEMPERATURE:
  //double P,T;
  
  void auto_0(){

    //CMAX = env::CMAX;
    //OMAX = env::OMAX;

    //P = env::P;
    //T = env::T;
  }

  // FIND CARBON AND OXYGEN NUMBERS OF CELL:
  //=================================================
  // CARBON AND
  // OXYGEN NUMBERS:
  int CNO,ONO;
  
  void auto_1(){
    tie(CNO,ONO) = calc_cno_ono(ID,CMAX,OMAX);
  };
    
  // CALCULATE MOLECULAR WEIGHT:
  //=================================================
  // MOLECULAR WEIGHT:
  double MW;
  
  void auto_2(){
    MW = double(CNO)*12.0107 + double(ONO)*15.999 + ((double(CNO)*2) + 2 - double(ONO))*1.00794;
  };

  // CALCULATE REACTIVITY WITH OH:
  //=================================================
  // kOH [cm3 molecules-1 s-1] AND
  // kOH IN PPM [ppm-1 s-1]:
  double kOH,kOH_ppm;
  
  void auto_3(){
    // BASED ON PARAMETERIZATIONS
    kOH = calc_koh(CNO,ONO,T);

    if ((CNO == 11)&&(ONO == 6)){
      kOH = 0.0;
    }

    if ((CNO == 11)&&(ONO == 4)){
      kOH = 0.0;
    }    
    
    // AVOGADRO'S NUMBER [molecules mol-1] AND
    // GAS CONSTANT [J mol-1 K-1]:
    double Na = 6.02e23;
    double R  = 8.314;

    // CONVERT TO PPM: 
    kOH_ppm = kOH*(7.3395e15)/300.;

    // 1ST-GEN. ONLY:
    if (G1ONLY == 1){
      kOH_ppm = 0.0;
    }
  };

  // THIS METHOD ALLOWS KOH TO BE CHANGED:
  //=================================================
  void change_koh(double kOH_ppm_in){
    kOH_ppm = kOH_ppm_in;
  }
  
  // FIND NUMBER OF FUNC. AND FRAG. DEPENDANTS:
  //=================================================
  // NUMBER OF FUNC. AND FRAG. DEPENDANTS:
  int n_func,n_frag;
  
  void auto_4(){
    n_func = find_nfunc(ID,CNO,ONO,CMAX,OMAX);
    n_frag = find_nfrag(ID,CNO,ONO,CMAX,OMAX);

    //if ((CNO == 10)&&(ONO == 7)){
    //  n_frag = n_frag + 1;
    //}
    
  };
  
  // FIND FUNC. DEPENDANTS:
  //=================================================
  // FUNC. DEPENDANTS:
  vector<int> depend_func;
  
  void auto_5(){
    depend_func = find_depend_func(ID,n_func);
  };
  
  // FIND FRAG. DEPENDANTS:
  //=================================================
  // FRAG. DEPENDANTS:
  vector<int> depend_frag;
  
  void auto_6(){
    depend_frag = find_depend_frag(ID,CNO,ONO,CMAX,OMAX,n_frag);

    //if ((CNO == 10)&&(ONO == 7)){
      //depend_frag.push_back(67);
      //depend_frag.push_back(66);
    //}
  };
  
  // EXTRACT THE PARAMETERS:
  //=================================================
  // PROBABILITIES OF FRAGMENTATION AND
  // FUNCTIONALIZATION:
  double pfrag,pfunc;

  // OXYGEN-ADDITION PROBABILITIES AND
  // THE SUM:
  vector<double> pfs;
  vector<double> pfsums;
  
  void auto_7(){

    double mfrag = params[0];
    pfrag = min(pow(double(ONO)/double(CNO),mfrag),1.);
    pfunc = 1. - pfrag;

    pfs.push_back(params[1]);
    pfs.push_back(params[2]);
    pfs.push_back(params[3]);
    pfs.push_back(params[4]);
  }
  
  // READ GAS-PHASE CONCENTRATION AT START TIME:
  //=================================================
  // INITIAL GAS CONCENTRATION [ppm]:
  double Cgas = 0.0;
  
  void read_cgas(double Cgas_in){
    Cgas = Cgas_in;
  }
  
  // INITIALIZE ADD AND DEL VARIABLES:
  //=================================================
  double to_add = 0.;
  double to_del = 0.;

  // THIS METHOD STEPS FOR GAS-PHASE CHEMISTRY:
  //=================================================
  void step_chem(double OH_in,double dt,vector<SOM_CELLS> &list_cells){

    // RADICAL CONCENTRATION:
    double OH;
    OH = OH_in;
    
    // CALCUATE REACTED AMOUNT (RUNGE-KUTTA):
    double k1,k2,k3,k4;

    k1 = kOH_ppm*OH*Cgas*(dt);
    k2 = kOH_ppm*OH*(Cgas+0.5*k1)*(dt);
    k3 = kOH_ppm*OH*(Cgas+0.5*k2)*(dt);
    k4 = kOH_ppm*OH*(Cgas+k3)*(dt);
    
    to_del = to_del + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
    
    // ADD TO FUNC. DEPENDANTS:
    int ID_add;
    
    for (int i = 0;i <= n_func-1;i++){

      ID_add = depend_func[i];

      // GET THE CELL:
      SOM_CELLS &cell_add = list_cells[ID_add-1];

      // NOW ADD TO CELL:
      double pf;{
	if (i == n_func-1){
	  pf = accumulate(pfs.begin()+i,pfs.end(),0.);
	} else {
	  pf = pfs[i];
	}
      }
      cell_add.to_add = cell_add.to_add + to_del*pfunc*pf;
    }
    
    // KEEP MASS IN CELL FOR UPPER BOUND CELLS:
    if (ONO == min(CNO*2,OMAX)){      
      to_add += to_del*pfunc;
    }
    
    // ADD TO FRAG. DEPENDANTS:
    for (int i = 0;i <= n_frag-1;i++){

      ID_add = depend_frag[i];

      // GET THE CELL:
      SOM_CELLS &cell_add = list_cells[ID_add-1];
      
      // NOW ADD TO CELL:
      double fac = 2./double(n_frag);
      fac = roundf(fac*1000)/1000;
      
      cell_add.to_add = cell_add.to_add + to_del*pfrag*fac;
    }
 
  };

  // SUM UP ADDITION AND REDUCTION:
  //=================================================
  void sum_chem(){
    Cgas   = Cgas + to_add - to_del;
    to_add = 0.;
    to_del = 0.;
  }
  

};

