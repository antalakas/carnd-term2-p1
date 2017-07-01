#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * A helper method to calculate polar coordinates from cartesian.
  */
  VectorXd CartesianToPolar(const VectorXd& x_state);

  /**
  * A helper method to calculate cartesian coordinates from polar.
  */
  VectorXd PolarToCartesian(const float& ro, const float& phi);

  /**
  * A helper method to normalize phi in the range -Pi -> Pi.
  * The Kalman Filter expects low values for phi.
  */
  float NormalizePhi(const float& phi);
};

#endif /* TOOLS_H_ */
