//==================================================================
// THIS FUNCTION FINDS THE NUMBER OF DEPENDANTS OF FUNCTIONALIZATION
//==================================================================

#pragma once

int find_nfunc(int ID,int CNO,int ONO,int CMAX,int OMAX){

  // FIND THE NUMBER OF FUNC. DEPENDANTS:
  //================================================================  
  // THE UPPER BOUND OF FUNC. PRODUTS:
  int cap;{
    cap = min({CNO*2,OMAX,ONO+4});
  }

  // FIND THE NUMBER OF FUNC. PRODUCTS:
  int n_func;{
    n_func = cap - ONO;
  }
  
  return n_func;
}
