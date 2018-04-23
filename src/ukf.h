#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "tools.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;


class UKF {
public:

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * handles everything after a new incoming measurement, including
   * prediction and update
   * @param {MeasurementPacka
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);


private:

  Tools tools;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state covariance matrix
  MatrixXd P_;

  /// * measurement covariance matrix
  MatrixXd R_;

  ///* measurement covariance matrix for laser
  MatrixXd R_laser_;

  ///* measurement covariance matrix for radar
  MatrixXd R_radar_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long previous_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_phidd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  /// * Number of sigma points
  int n_sig_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* measurement dimension for laser
  int n_z_laser_;

  ///* measurement dimension for radar
  int n_z_radar_;

  ///* predicted measurement
  VectorXd z_pred_;

  ///* predicted sigma points in measurement space
  MatrixXd Zsig_;

  ///* predicted measurement covariance
  MatrixXd S_;

  /**
   * Initialize with the first measurement
   * @param {MeasurementPackage} meas_package The latest measurement data of
   * either radar or laser.
   */
  void Initialize(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Generates sigma points
   * Returns matrix with number of cols as number of sigma points
   * and each row represents a sigma point
   */
  MatrixXd GenerateSigmaPoints();

  /**
   * Predicts sigma points as a matrix of n_x_ * n_sig_,
   * with number of cols as number of sigma points
   * and each row represents a predicted sigma point with dimension n_x_
   */
  void PredictSigmaPoints(MatrixXd Xsig, double delta_t);

  /**
   * Predicts state mean and covariance using weighted sum
   * Updates state vector x_ and state covariance matrix P_
   */
  void PredictState(MatrixXd X_pred_);

  /**
   * Wrapper function that updates the state and the state covariance matrix
   * @param meas_package The measurement at k+1
   */
  void Update(MeasurementPackage meas_package);

  /* Convert predicted state mean and covariance to 
   * a measurement mean and covariance
   * Result: updates S_ and z_pred
   */
  void UpdateMeasurement(int sensor_type_);

  /**
  * Updates the state and the state covariance matrix using 
  * predicted measurement mean and covariance of the previous measurement
  * as well as the latest raw measurement
  * @param {VectorXd} z: latest raw measurement
  */
  void UpdateState(VectorXd z);
};

#endif /* UKF_H */
