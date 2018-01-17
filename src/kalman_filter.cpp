#include "kalman_filter.h"
#include <iostream>


using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    DONE: * predict the state
  */
  
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::NormalizeAngle(double& phi){
  /**
    * Using the atan2 function, normalize the angle phi so that it lies
    * between (-pi/2,pi/2)
   */

  phi = atan2(sin(phi), cos(phi));
}

void KalmanFilter::KF(const VectorXd& y){
  /**
   * Perform the Linear Algebra calculations, common for the Update() and
   * the UpdateEKF() functions
   */

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
    
  //new estimate
  x_ = x_ + K * y;
  P_ -= K * H_ * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    DONE: * update the state by using Kalman Filter equations
  */
  
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  
  //Call Estimator function to complete the update
  KF(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    DONE: * update the state by using Extended Kalman Filter equations
  */
  
  //recover state parameters
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  
  //Precalculate reused quantities for error calculation
  double c1 = px*px + py*py;
  double c2 = sqrt(c1);
  
  //check division by zero
  if(fabs(c1) < 0.0001){
    std::cout << "UpdateEKF () - Error - Divide by zero,  skipping current RADAR update calculation.\n\n";
    return;
  }
  
  VectorXd h_x = VectorXd(3);
  h_x << c2, atan2(py, px), (px*vx + py*vy)/c2;
 
  VectorXd y = z - h_x;
  std::cout << std::endl << std::endl <<"Phi = " << y(1) << std::endl;

  //Adjust phi value to be between -pi and pi.
  NormalizeAngle(y(1));
  
  //Call Estimator function to complete the update
  KF(y);
  
}
