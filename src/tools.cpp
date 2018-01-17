#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    DONE: * Calculate the RMSE here.
  */
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  int n = estimations.size();
  int g = ground_truth.size();
  //cout << "number of estimations = " << n << endl;
  //cout << "number of ground truths = " << g << endl;
  
  if(n==0 || n!=g) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }
  
  //accumulate squared residuals
  for(int i=0; i < n; ++i) {
    
      //cout << "i = " << i << endl;
      VectorXd res = estimations[i]-ground_truth[i];
      res = res.array() * res.array();
      rmse += res;
      //cout << "residuals = " << residuals << endl << "sum = " << rmse << endl;
  }
  
  //calculate the mean
  rmse = rmse/n;
  //cout << "mean = " << rmse << endl;
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  //cout << "RMSE = " << rmse << endl;
  
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    DONE: * Calculate a Jacobian here.
  */
  
  MatrixXd Hj(3,4);
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  //Precalculate parameters
  double c1 = px*px + py*py;
  double c2 = sqrt(c1);
  double c3 = c1*c2;
  
  //check for division by zero
  c1 = std::max(0.001, c1);
 
  //compute the Jacobian matrix
  Hj << (px/c2),                 (py/c2),                 0,        0,
        -(py/c1),                (px/c1),                 0,        0,
        py*(vx*py - vy*px)/c3,    px*(vy*px - vx*py)/c3,  px/c2,    py/c2;
  
  return Hj;
}
