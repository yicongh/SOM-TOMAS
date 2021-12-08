//==============================================================
// THIS FUNCTION FINDS THE NUMBER OF DEPENDANTS OF FRAGMENTATION
//==============================================================

#pragma once

int find_nfrag(int ID,int CNO,int ONO,int CMAX,int OMAX){
  
  // FIND THE NUMBER OF FRAG. DEPENDANTS:
  //============================================================
  // BASE NUMBER:
  int BASE_N;{
    BASE_N = (CNO - 1)*(ONO + 1);
  }
  
  // NUMBER OF POSITIONS OUT OF LIMITS:
  int OUT_N = 0;

  for (int i = 1; i <= CNO-1; i++){
    for (int j = 1; j <= ONO+1; j++){
      if (j > 2*i){
	OUT_N = OUT_N + 2;
      } else if (j > OMAX){
	OUT_N = OUT_N + 1;
      }
    }
  }

  // NUMBER OF FRAG. DEPENDANTS:
  int n_frag;{
    n_frag = BASE_N - OUT_N;

    if (n_frag < 0){
      n_frag = 0;
    }
  }
   
  return n_frag;

}
