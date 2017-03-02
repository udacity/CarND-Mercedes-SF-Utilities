/**
 * ukf.cpp
 * @author MBRDNA
 *
 */

#include <iostream>
#include "ukf.h"

UKF::UKF() {
  // TODO Auto-generated constructor stub
  Init();
}

UKF::~UKF() {
  // TODO Auto-generated destructor stub
}

void UKF::Init() {
  std::cout << "Init" << std::endl;

  x_ = VectorXd(5);

  // initial position
  x_ << 0.1, 0.1, 0.1, 0.1, 0.01;

  P_ = MatrixXd(5, 5);

  //initial state covariance matrix
  P_ << 0.2, 0, 0, 0, 0,
        0, 0.2, 0, 0, 0,
        0, 0, 0.2, 0, 0,
        0, 0, 0, 0.3, 0,
        0, 0, 0, 0, 0.3;

  time_us_ = 1;

  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = kappa_/(kappa_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  // 2n+1 weights
    double weight = 0.5/(n_aug_+kappa_);
    weights_(i) = weight;
  }

  Xi_pred_ = MatrixXd(n_, 2 * n_aug_ + 1);
  Xi_pred_.fill(0.0);
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
    }
    return;
  }

  // check if the incoming sensor type is to be used
  if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
      return;

  if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
      return;

  // calculate the time difference of the measurement to the current state
  long delta_t_us = meas_package.timestamp_ - time_us_;
  double delta_t= static_cast<double>(delta_t_us) / static_cast<double>(1e6) ;
  time_us_ = meas_package.timestamp_;

  std::cout << "delta_t = " << delta_t << std::endl;
  std::cout << "meas_package.timestamp_ = " << meas_package.timestamp_ << std::endl;
  std::cout << "meas_package.sensor_type_ = " << meas_package.sensor_type_ << std::endl;


  // call the prediction function
  Prediction(delta_t);

  // distinguish between radar measurement and laser measurement for update call
  if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else
  {

  }

}

void UKF::Prediction(double delta_t) {
  std::cout << "Prediction" << std::endl;

  // build augmented vector
  VectorXd x_aug = VectorXd(7);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // build augmented state covariance
  MatrixXd P_aug_ = MatrixXd(7, 7);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  MatrixXd L = P_aug_.llt().matrixL();

  //std::cout << "P_aug_ = " << std::endl << P_aug_ << std::endl;
  std::cout << "L = " << std::endl << L << std::endl;

  // create sigma point matrix
  MatrixXd Xi_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  double sqrt_npk = sqrt(n_aug_+kappa_);
  Xi_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xi_aug.col(i+1)        = x_aug + sqrt_npk * L.col(i);
    Xi_aug.col(i+1+n_aug_)  = x_aug - sqrt_npk * L.col(i);
  }

  //std::cout << "Xi_aug = " << std::endl << Xi_aug << std::endl;

  // predict sigma points
  Xi_pred_ = MatrixXd(n_, 2 * n_aug_ + 1);
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    // extract values for better readability
    double p1 = Xi_aug(0,i);
    double p2 = Xi_aug(1,i);
    double v = Xi_aug(2,i);
    double yaw = Xi_aug(3,i);
    double yawd = Xi_aug(4,i);
    double nu_a = Xi_aug(5,i);
    double nu_yawdd = Xi_aug(6,i);

    double p1_p, p2_p;

    if (abs(yawd) > 0.001) {
        p1_p = p1 + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        p2_p = p2 + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        p1_p = p1 + v*delta_t*cos(yaw);
        p2_p = p2 + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    p1_p = p1_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    p2_p = p2_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xi_pred_(0,i) = p1_p;
    Xi_pred_(1,i) = p2_p;
    Xi_pred_(2,i) = v_p;
    Xi_pred_(3,i) = yaw_p;
    Xi_pred_(4,i) = yawd_p;
  }

  std::cout << "Xi_pred_ = " << std::endl << Xi_pred_ << std::endl;

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    x_ = x_+ weights_(i) * Xi_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    P_ = P_ + weights_(i) * ( Xi_pred_.col(i) - x_ ) * ( Xi_pred_.col(i) - x_ ).transpose() ;
  }

  std::cout << "Predicted state ukf" << std::endl;
  std::cout << x_ << std::endl;

  std::cout << "Predicted covariance matrix ukf" << std::endl;
  std::cout << P_ << std::endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  std::cout << "UpdateLidar" << std::endl;

  MatrixXd Xi_z = Xi_pred_.topRows(2); // laser can measure p1 and p2

  std::cout << "Xi_z: " << Xi_z << std::endl;


  // predict measurement
  VectorXd z_pred = VectorXd(2);
  z_pred.fill(0.0);

  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Xi_z.col(i);
  }
  std::cout << "z_pred: " << z_pred << std::endl;


  // measurement covariance matrix S
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    S = S + weights_(i) * ( Xi_z.col(i) - z_pred ) * ( Xi_z.col(i) - z_pred ).transpose();
  }
  // add measurement noise
  MatrixXd R = MatrixXd(2,2);
  R << std_lasp1_*std_lasp1_,                     0,
                           0, std_lasp2_*std_lasp2_;
  S = S + R;

  std::cout << "S: " << S << std::endl;
  std::cout << "weights_: " << weights_ << std::endl;


  // cross correlation Tc
  MatrixXd Tc = MatrixXd(n_,2);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    Tc = Tc + weights_(i) * ( Xi_pred_.col(i) - x_ ) * ( Xi_z.col(i) - z_pred ).transpose();
  }

  std::cout << "Tc: " << Tc << std::endl;


  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  std::cout << "K: " << K << std::endl;

  // update
  VectorXd z = meas_package.raw_measurements_;
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K*S*K.transpose();

  std::cout << "z: " << z << std::endl;

  std::cout << "Updated state x_: " << x_ << std::endl;
  std::cout << "Updated state covariance P_: " << std::endl;
  std::cout <<  P_ << std::endl;

/*  //Error y
  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - H_ * x_;

  //Kalman Gain
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;

  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //update
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(2, 2);
  P_ = (I - K * H_) * P_;*/
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {


  std::cout << "UpdateRadar" << std::endl;

  // transform sigma points into measurement space
  MatrixXd Xi_z = MatrixXd(3, 2 * n_aug_ + 1); // radar can measure r, phi, and r_dot
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    double p_x = Xi_pred_(0,i);
    double p_y = Xi_pred_(1,i);
    double v  = Xi_pred_(2,i);
    double yaw = Xi_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    Xi_z(0,i) = sqrt(p_x*p_x + p_y*p_y);                            //r
    Xi_z(1,i) = atan2(p_y,p_x);                                 //phi
    Xi_z(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);       //r_dot
  }


  std::cout << "Xi_z: " << Xi_z << std::endl;

  // predict measurement
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0.0);

  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Xi_z.col(i);
  }
  std::cout << "z_pred: " << z_pred << std::endl;


  // measurement covariance matrix S
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    S = S + weights_(i) * ( Xi_z.col(i) - z_pred ) * ( Xi_z.col(i) - z_pred ).transpose();
  }
  // add measurement noise
  MatrixXd R = MatrixXd(3,3);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  std::cout << "S: " << S << std::endl;
  std::cout << "weights_: " << weights_ << std::endl;


  // cross correlation Tc
  MatrixXd Tc = MatrixXd(n_,3);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    Tc = Tc + weights_(i) * ( Xi_pred_.col(i) - x_ ) * ( Xi_z.col(i) - z_pred ).transpose();
  }

  std::cout << "Tc: " << Tc << std::endl;


  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  std::cout << "K: " << K << std::endl;

  // update
  VectorXd z = meas_package.raw_measurements_;
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K*S*K.transpose();

  std::cout << "z: " << z << std::endl;

  std::cout << "Updated state x_: " << x_ << std::endl;
  std::cout << "Updated state covariance P_: " << std::endl;
  std::cout <<  P_ << std::endl;




}


