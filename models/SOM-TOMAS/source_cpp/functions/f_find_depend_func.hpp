//========================================================
// THIS FUNCTION FINDS THE DEPENDANTS OF FUNCTIONALIZATION
//========================================================

#pragma once

vector<int> find_depend_func(int ID,int n_func){
  
  // FIND FUNC. DEPENDANT IDS:
  //======================================================
  // LIST OF DEPENDANT IDS:
  vector<int> depend_func;{
    for (int i = 1; i <= n_func; i++){
      depend_func.push_back(ID + i);
    }
  }
  return depend_func;
}
