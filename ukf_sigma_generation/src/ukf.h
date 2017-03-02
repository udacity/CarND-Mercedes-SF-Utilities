/**
 * ukf.h
 * @author MBRDNA
 */

#ifndef UKF_H
#define UKF_H
#include "../../Eigen/Dense"
#include "measurement_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_ = false;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_ = true;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_ = true;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	VectorXd x_;

	///* state covariance matrix
	MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xi_pred_;

  ///* time when the state is true, in us
  long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_ = 0.2;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_ = 0.2;

  ///* Laser measurement noise standard deviation position1 in m
  double std_lasp1_ = 0.15;

  ///* Laser measurement noise standard deviation position2 in m
  double std_lasp2_ = 0.15;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_ = 0.3;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_ = 0.0175;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ = 0.1;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_ = 5;

  ///* Augmented state dimension
  int n_aug_ = 7;

  ///* Sigma point spreading parameter
  double kappa_ = 3-n_aug_;

  /**
	 * Constructor
	 */
	UKF();

	/**
	 * Destructor
	 */
	virtual ~UKF();

	/**
	 * Init Initializes Unscented Kalman filter
	 */
	void Init();

	/**
	 * ProcessMeasurement
	 */
	void ProcessMeasurement(MeasurementPackage meas_package);

	/**
	 * Prediction Predicts sigma points, the state, and the state covariance
	 * matrix
	 * @param delta_t Time between k and k+1 in s
	 */
	void Prediction(double delta_t);

	/**
	 * Updates the state and the state covariance matrix using a laser measurement
	 * @param meas_package The measurement at k+1
	 */
	void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);




};

#endif /* UKF_H */


