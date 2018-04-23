#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_phidd_ = 1.;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Weights of sigma points
  weights_ = VectorXd(n_sig_);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial state covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);
  P_(0, 0) = 1.;
  P_(1, 1) = 1.;
  P_(2, 2) = 1.;
  P_(3, 3) = 1.;
  P_(4, 4) = 1.;

  ///* measurement dimension for laser
  n_z_laser_ = 2;

  ///* measurement dimension for radar
  n_z_radar_ = 3;

  // measurement covariance matrix for laser
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_.fill(0.0);
  R_laser_(0, 0) = pow(std_laspx_, 2);
  R_laser_(1, 1) = pow(std_laspy_, 2);

  // measurement covariance matrix for radar
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_.fill(0.0);
  R_radar_(0, 0) = pow(std_radr_, 2);
  R_radar_(1, 1) = pow(std_radphi_, 2);
  R_radar_(2, 2) = pow(std_radrd_, 2);
}

UKF::~UKF() {}

/**
 * handles everything after a new incoming measurement, including
 * prediction and update
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_){
    Initialize(meas_package);
    return;
  }
  
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER &&
      use_laser_) ||
        (meas_package.sensor_type_ == MeasurementPackage::RADAR &&
      use_radar_)){

    double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;;
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(delta_t);
    Update(meas_package);
  }
}

/**
 * Initialize with the first measurement
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::Initialize(MeasurementPackage meas_package){
  // if radar, assume phid is 0
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    float r = meas_package.raw_measurements_[0];
    float phi = meas_package.raw_measurements_[1];
    float rd = meas_package.raw_measurements_[2];
    x_ << r * cos(phi), r * sin(phi), rd, phi, 0;
  }
  // if laser, assume v, phi and phid are all 0
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
  }
  previous_timestamp_ = meas_package.timestamp_;
  is_initialized_ = true;
  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_ = GenerateSigmaPoints();
  PredictSigmaPoints(Xsig_, delta_t);
  PredictState(Xsig_pred_);
}

/**
 * Generates sigma points
 * Returns matrix of n_aug_ * n_sig_,
 * with number of cols as number of sigma points
 * and each row represents a sigma point
 */
MatrixXd UKF::GenerateSigmaPoints(){
  //create sigma point matrix
  MatrixXd Xsig_ = MatrixXd(n_aug_, n_sig_);

  //create augmented mean vector
  VectorXd x_aug_ = VectorXd(7);
  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_, n_x_) = std_a_ * std_a_;
  P_aug_(n_x_ + 1, n_x_ + 1) = std_phidd_ * std_phidd_;

  //create square root matrix
  MatrixXd A = P_aug_.llt().matrixL();
  
  //create augmented sigma points
  Xsig_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++){
      Xsig_.col(1 + i) = x_aug_ + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_.col(1 + n_aug_ + i) = x_aug_ - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  return Xsig_;
}

/**
 * Predicts sigma points, as a matrix of n_x_ * n_sig_,
 * with number of cols as number of sigma points
 * and each row represents a predicted sigma point with dimension n_x_
 * Updates Xsig_pred_
 */
void UKF::PredictSigmaPoints(MatrixXd Xsig_, double delta_t){
  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  Xsig_pred_.fill(0.0);

  for (int i = 0; i < n_sig_; i++){
      VectorXd xk_aug = Xsig_.col(i); 
      float px = xk_aug(0);
      float py = xk_aug(1);
      float v = xk_aug(2);
      float phi = xk_aug(3);
      float phi_dot = xk_aug(4);
      float nu_a = xk_aug(5);
      float nu_phidd = xk_aug(6);

      VectorXd xk = VectorXd(5);
      xk << px, py, v, phi, phi_dot;

      VectorXd integral = VectorXd(5);
      integral << v * cos(phi) * delta_t,
                  v * sin(phi) * delta_t,
                  0,
                  phi_dot * delta_t,
                  0;
      if (fabs(phi_dot) > EPSILON){
          integral[0] = v / phi_dot * (sin(phi + phi_dot * delta_t) - sin(phi));
          integral[1] =  v / phi_dot * (-cos(phi + phi_dot * delta_t) + cos(phi));
      }

      VectorXd noise = VectorXd(5);
      noise << 1./2 * delta_t * delta_t * cos(phi) * nu_a,
               1./2 * delta_t * delta_t * sin(phi) * nu_a,
               delta_t * nu_a,
               1./2 * delta_t * delta_t * nu_phidd,
               delta_t * nu_phidd;
      
      Xsig_pred_.col(i) = xk + integral + noise;
  }
}

/**
 * Predicts state mean and covariance using weighted sum
 * Updates state vector x_ and state covariance matrix P_
 */
void UKF::PredictState(MatrixXd Xsig_pred_){
  x_ = tools.GetWeightedMean(weights_, Xsig_pred_);
  P_ = tools.GetWeightedCovariance(weights_, Xsig_pred_, x_, 3);
}

/**
 * Wrapper function that updates the state and the state covariance matrix
 * incorporating the latest measurement
 * @param meas_package The measurement at k+1
 */
void UKF::Update(MeasurementPackage meas_package){
  UpdateMeasurement(meas_package.sensor_type_);
  UpdateState(meas_package.raw_measurements_);
}

/**
 * Convert predicted state mean and covariance to 
 * a measurement mean and covariance
 * Result: updates S_ and z_pred_
 */
void UKF::UpdateMeasurement(int sensor_type_){
  ArrayXd px = Xsig_pred_.row(0);
  ArrayXd py = Xsig_pred_.row(1);
  ArrayXd v = Xsig_pred_.row(2);
  ArrayXd phi = Xsig_pred_.row(3);
  ArrayXd phi_dot = Xsig_pred_.row(4);

  if (sensor_type_ == MeasurementPackage::LASER){
    Zsig_ = MatrixXd(n_z_laser_, n_sig_);
    Zsig_.fill(0.0);
    Zsig_.row(0) = px;
    Zsig_.row(1) = py;
  }
  else{
    Zsig_ = MatrixXd(n_z_radar_, n_sig_);
    Zsig_.fill(0.0);
    Zsig_.row(0) = (px * px + py * py).cwiseSqrt();
    for (int i = 0; i < n_sig_; i ++){
      Zsig_(1, i) = atan2(py(i), px(i));
    }
    ArrayXd denom = Zsig_.row(0);
    Zsig_.row(2) = (px * cos(phi) + py * sin(phi)) * v / denom;
  }
  
  //calculate mean predicted measurement
  z_pred_ = tools.GetWeightedMean(weights_, Zsig_);

  //calculate innovation covariance matrix S
  S_ = MatrixXd(n_z_laser_, n_z_laser_);
  S_.fill(0.0);

  int angle_index =  (sensor_type_ == MeasurementPackage::RADAR) ? 1 : -1;
  S_ = tools.GetWeightedCovariance(weights_, Zsig_, z_pred_, angle_index);
  
  if (sensor_type_ == MeasurementPackage::LASER){
    S_(0, 0) += pow(std_laspx_, 2);
    S_(1, 1) += pow(std_laspy_, 2);
  }
  else{
    S_(0, 0) += pow(std_radr_, 2);
    S_(1, 1) += pow(std_radphi_, 2);
    S_(2, 2) += pow(std_radrd_, 2);
  }
}

/**
 * Updates the state and the state covariance matrix using 
 * predicted measurement mean and covariance of the previous measurement
 * as well as the latest raw measurement
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateState(VectorXd z) {
  MatrixXd T = weights_[0] * (Xsig_pred_.col(0) - x_
      ) * (Zsig_.col(0) - z_pred_).transpose();
  for (int i = 1; i < n_sig_; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.NormalizeAngle(x_diff(3));
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    z_diff(1) = tools.NormalizeAngle(z_diff(1));
    T += weights_[i] * x_diff * z_diff.transpose();
  }
  MatrixXd K = T * S_.inverse();
  x_ += K * (z - z_pred_);
  P_ -= K * S_ * K.transpose();
}
