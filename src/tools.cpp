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
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0){
      cout << "Invalid estimation or ground_truth data" << endl;
      return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
      VectorXd res = estimations[i] - ground_truth[i];

      res = res.array() * res.array();

      rmse += res;
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by zero
  if(fabs(c1) < 0.0001){
      cout << "CalculateJacobian () - Error - Division by Zero" << endl;
      return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}

VectorXd Tools::CartesianToPolar(const VectorXd &x_state) {
  /**
  TODO:
    * Calculate polar coordinates here.
  */
  VectorXd hx(3);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px*px+py*py;
//  float c2 = py/px;
  float c3 = px*vx + py*vy;

  //check division by zero
  if(fabs(c1) < 0.0001){
      cout << "CartesianToPolar () - Error - Division by Zero" << endl;
      return hx;
  }

  float h1 = sqrt(c1);
  float h2 = atan2(py, px);
  float h3 = c3/c1;

  hx << h1, h2, h3;

  return hx;
}

/**
* A helper method to calculate cartesian coordinates from polar.
*/
VectorXd Tools::PolarToCartesian(const float& ro, const float& phi)
{
  /**
  TODO:
    * Calculate cartesian coordinates here.
  */
  VectorXd coords(2);

  coords << ro * cos(phi), ro * sin(phi);

//  float tan2 = tan(phi) * tan(phi);
//
//  float c1 = 1 + tan2;
//
//  //check division by zero
//  if(fabs(c1) < 0.0001){
//      cout << "PolarToCartesian () - Error - Division by Zero" << endl;
//      return coords;
//  }
//
//  float c2 = 1/c1;
//
//  float px = ro * sqrt(c2);
//  float py = ro * sqrt(1-c2);
//
//  coords << px, py;

  return coords;
}

float Tools::NormalizePhi(const float& phi)
{
  /**
  TODO:
    * Normalize phi here.
  */
  return fmod(phi, M_PI);
}