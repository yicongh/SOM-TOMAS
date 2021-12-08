//====================================================
// THIS FUNCTION FINDS THE DEPENDANTS OF FRAGMENTATION
//====================================================

#pragma once

vector<int> find_depend_frag(int ID,int CNO,int ONO,int CMAX,int OMAX,int n_frag){
  
  // FIND FRAG. DEPENDANT IDS:
  //==================================================
  // LIST OF FRAG. DEPENDANT IDS:
  vector<int> depend_frag;
  
  // COMPLEMENTARY POSITIONS:
  int ix,jx;
  
  // LOOP OVER THE BASE GRID:
  for (int i = 1; i <= CNO-1; i++){
    for (int j = 1; j <= ONO+1; j++){

      // COUNTER:
      int ctr = 0;
      
      // CALCULATE COMPLEMENTARY POSITIONS:
      ix = (CNO-1) + 1 - i;
      jx = (ONO+1) + 1 - j;

      // FOR OUT OF LIMIT POSITIONS:
      if (j > 2*i || j > OMAX){
	continue;
      } else {
	for (int k = 1; k <= i-1; k++){
	  ctr = ctr + min(2*k,OMAX);
	}
	ctr = ctr + j;
      }
      
      // SKIP THE ABSENT POSITIONS:
      if (jx > 2*ix);
      else{
	depend_frag.push_back(ctr);
      }
    }
  }
  return depend_frag;
}
