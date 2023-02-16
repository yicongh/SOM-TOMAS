//==========================================================
// This function calculates the cell ID based on CNO and ONO
//==========================================================

#pragma once

int calc_id(int CNO,int ONO,int CMAX,int OMAX){

  // Stack height of each grid column:
  int height[CMAX];

  for (int i = 0; i <= CMAX - 1; i = i + 1){

    height[i] = min(OMAX,2*(i+1));
  }

  // Cumulative height at each grid column:
  int cumu_h[CMAX];

  partial_sum(height,height+CMAX,cumu_h);

  // Find the ID of the cell:
  int ID;{

    // First part of ID due to CNO:
    int ID_p1;
    ID_p1 = accumulate(height,height+CNO-1,0);

    // Second part of ID due to ONO:
    int ID_p2;
    ID_p2 = ONO;

    // Combine the two:
    ID = ID_p1 + ID_p2;
  }
  return ID;
}
