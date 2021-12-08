//=====================================================
// This function calculates the kOH of the grid species
//=====================================================

#pragma once

double calc_koh(int CNO,int ONO,double T){

  double A1,A2,A3;
  A1 = -15.103;
  A2 = -3.9481;
  A3 = -0.79796;

  double b1,b2;
  b1 = -0.2583*double(CNO) + 5.8944;

  if (CNO <= 15){
    b2 = 0.0314*double(CNO) + 0.9871;} else {
    b2 = 0.2500*double(CNO) - 2.1830;
  }

  double sigma;
  if (CNO <= 15){
    sigma = 0.0214*double(CNO) + 0.5238;} else{
    sigma = -0.115*double(CNO) + 2.695;
  }

  double Ea = 1000.0;

  double kOH;{
    
    double p1,p2,p3;

    p1 = pow(10.,A1 + A2 * pow(double(CNO),A3));

    p2 = 1.0 + b1/sigma/sqrt(2.0*3.14159265) * \
         exp(-1.0*pow(log(double(ONO)+0.01)-log(b2),2.)/2./pow(sigma,2.));

    p3 = pow(298.,2.) * exp(-1.0*Ea/8.314/298.);

    kOH = p1*p2*p3;
  }

  // cout << CNO << ONO << "kOH = " << kOH << endl;
  
  return kOH;

}
