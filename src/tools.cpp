#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse.fill(0.0);

  if (estimations.size() != ground_truth.size()){
    cout << "estimations is of size " << estimations.size()
         << " but ground_truth is of size " << estimations.size() << endl;
    return rmse;
  }
  if (estimations.size() == 0){
    cout << "estimations has size 0." << endl;
    return rmse;
  }

  for (int i = 0; i < estimations.size(); i++){
    VectorXd e = estimations[i];
    VectorXd g = ground_truth[i];
    VectorXd squared = (e - g).array() * (e - g).array();
    rmse += squared;
  }

  rmse /= estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

/**
 * + / - 2 PI from an angle until it is between -PI and PI
 * @param {double} yaw: the angle
 */ 
double Tools::NormalizeAngle(double phi){
  if (phi > 0) {
    return fmod(phi + PI, 2 * PI) - PI;
  } else {
    return fmod(phi - PI, 2 * PI) + PI;
  }
}

/**
 * get weighted mean as sum of weight * m_i
 * where m_i is a column in m
 * @param {VectorXd} weights, dimension is of number of cols in m
 * @param {MatrixXd} m: each column is a data point
 * Return: a vector representing the mean
 */
VectorXd Tools::GetWeightedMean(VectorXd weights, MatrixXd m){
    VectorXd x = VectorXd(m.rows());
    //predict state mean
    x.fill(0.);
    for (int i = 0; i < weights.size(); i++){
        x += weights[i] * m.col(i);
    }
    return x;
}

/**
 * get covariance as sum of weight * (m_i - mean) * (m_i - mean).T
 * where m_i is a column in m
 * @param {VectorXd} weights, dimension is of number of cols in m
 * @param {MatrixXd} m: each column is a data point
 * @param {VectorXd} mean: a vector representing the mean
 * @param {int} angle_index: index in the vector mean corresponding to
 *   an angle measurement, so we will normalize the angle
 *   input -1 if no angle measurement is involved
 * Return: a matrix representing the covariance
 */
MatrixXd Tools::GetWeightedCovariance(VectorXd weights, MatrixXd m,
    VectorXd mean, int angle_index){
  MatrixXd covariance = MatrixXd(mean.size(), mean.size());
  covariance.fill(0.);

  for (int i = 0; i < weights.size(); i++){
      VectorXd diff = m.col(i) - mean;
      if (angle_index != -1){
        diff(angle_index) = NormalizeAngle(diff(angle_index));
      }
      covariance += weights[i] * diff * diff.transpose();
  }
  return covariance;
}
