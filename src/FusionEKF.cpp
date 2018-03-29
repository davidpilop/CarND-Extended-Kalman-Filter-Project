#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define EPS 0.0001 // A very small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    double px, py, vx, vy;

    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "EKF : First measurement RADAR" << endl;
      // Convert radar from polar to cartesian coordinates and initialize state.
      double rho = measurement_pack.raw_measurements_[0]; // range
      double phi = measurement_pack.raw_measurements_[1]; // bearing
      double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho

      // Coordinates convertion from polar to cartesian
      px = rho * cos(phi);
      py = rho * sin(phi);
      vx = rho_dot * cos(phi);
      vy = rho_dot * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // No velocity and coordinates are cartesian already.
      cout << "EKF : First measurement LASER" << endl;
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
      vx = 0;
      vy = 0;
    }

    // Deal with the special case initialisation problems
    if ( px < EPS ) px = EPS;
    if ( py < EPS ) py = EPS;
    ekf_.x_ << px, py, vx, vy;

    // Print the initialization results
    cout << "EKF init: " << ekf_.x_ << endl;
    // Saving first timestamp in seconds
    previous_timestamp_ = measurement_pack.timestamp_ ;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
    * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
    * Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   // Calculate the timestep between measurements in seconds
   double dt = measurement_pack.timestamp_ - previous_timestamp_;
   dt /= 1000000.0; // convert micros to s
   previous_timestamp_ = measurement_pack.timestamp_;

  // State transition matrix update
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // Noise covariance matrix computation
  // Noise values from the task
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  double dt_2 = dt * dt; //dt^2
  double dt_3 = dt_2 * dt; //dt^3
  double dt_4 = dt_3 * dt; //dt^4
  double dt_4_4 = dt_4 / 4; //dt^4/4
  double dt_3_2 = dt_3 / 2; //dt^3/2
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
	           0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
	           dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
 	           0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
    * Use the sensor type to perform the update step.
    * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
  	ekf_.R_ = R_radar_;
  	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
  	ekf_.R_ = R_laser_;
  	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
