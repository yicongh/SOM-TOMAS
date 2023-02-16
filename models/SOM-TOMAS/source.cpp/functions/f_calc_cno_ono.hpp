//====================================================
// THIS FUNCTION CALCULATES CNO AND ONO FOR A GIVEN ID
//====================================================

#pragma once

tuple <int,int> calc_cno_ono(int ID,int CMAX,int OMAX){

  // FIND THE CARBON AND OXYGEN NUMBERS:
  //==================================================
  // CARBON AND
  // OXYGEN NUMBERS:
  int CNO,ONO;
  
  // COUNTER:
  int ctr = 0;
  
  // LOOP THROUGH GRID:
  for (int i = 1; i <= CMAX; i++){
    for (int j = 1; j <= OMAX; j++){

      // SKIP POSITIONS OUT OF LIMITS:
      if (j > 2*i){
	continue;
      } else {
	ctr = ctr + 1;
      }
      
      // EXIT WHEN COUNTER EQUALS ID:
      if (ctr == ID){
	CNO = i;
	ONO = j;
	break;
      }
    }
  }
  
  return make_tuple(CNO,ONO);
};
