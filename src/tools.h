#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

const float PI = 3.1415927;
const float EPSILON = 0.0001;

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
   * + / - 2 PI from an angle until it is between -PI and PI
   * @param {double} yaw: the angle
   */  
  double NormalizeAngle(double phi);

  /**
   * get weighted mean as sum of weight * m_i
   * where m_i is a column in m
   * @param {VectorXd} weights, dimension is of number of cols in m
   * @param {MatrixXd} m: each column is a data point
   * Return: a vector representing the mean
   */
  VectorXd GetWeightedMean(VectorXd weights, MatrixXd m);

  /**
   * get covariance as sum of weight * (m_i - mean) * (m_i - mean).T
   * where m_i is a column in m
   * @param {VectorXd} weights, dimension is of number of cols in m
   * @param {MatrixXd} m: each column is a data point
   * @param {VectorXd} mean: a vector representing the mean
   * @param {int} angle_index: index in the vector mean corresponding to
   *   an angle measurement, so we will normalize the angle
   * Return: a matrix representing the covariance
   */
  MatrixXd GetWeightedCovariance(VectorXd weights, MatrixXd m,
    VectorXd mean, int angle_index);
};

#endif /* TOOLS_H_ */