
#include <iostream>
#include "../../Eigen/Dense"
#include <vector>
#include "ukf.h"
#include "measurement_package.h"
#include "ground_truth_package.h"
#include <fstream>
#include <sstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

int main() {

	/*******************************************************************************
	 *  Local configuration			   	    										 *
	 *******************************************************************************/
	//bool use_gt = true;	//true = generate measurement packages together with GT
	/*******************************************************************************
	 *  Set Measurements															 *
	 *******************************************************************************/
	vector<MeasurementPackage> measurement_pack_list;
	vector<GroundTruthPackage> gt_pack_list;
	//
	//test list with measurement packages

	string in_file_name_ = "../matlab_examples/obj_pose-laser-radar-synthetic-ukf-input.txt";
	//string in_file_name_ = "../fusion_measurement_processor/data/obj_pose-laser-radar-ekf-input.txt";
	ifstream in_file(in_file_name_.c_str(), std::ifstream::in);

	string out_file_name_ = "data/obj_pose-laser-radar-ukf-output.txt";
	std::ofstream out_file_(out_file_name_.c_str(), std::ofstream::out);
	//out_file_.open(out_file_name_.c_str(),std::ofstream::out);
	if (!out_file_.is_open()) {
		cout << "Cannot Open input file: " << out_file_name_;
	}

	string line;

	while (getline(in_file, line)) {

		string sensor_type;
		MeasurementPackage meas_package;
		GroundTruthPackage gt_package;

		istringstream iss(line);
		iss >> sensor_type;	//reads first element from the current line
		long timestamp;
		if (sensor_type.compare("L") == 0) {	//laser measurement
			//read measurements
			meas_package.sensor_type_ = MeasurementPackage::LASER;
			meas_package.raw_measurements_ = VectorXd(2);
			float px;
			float py;
			iss >> px;
			iss >> py;
			meas_package.raw_measurements_ << px, py;
			iss >> timestamp;
			meas_package.timestamp_ = timestamp;
			measurement_pack_list.push_back(meas_package);

		} else if (sensor_type.compare("R") == 0) {
			//read measurements
			meas_package.sensor_type_ = MeasurementPackage::RADAR;
			meas_package.raw_measurements_ = VectorXd(3);
			float ro;
			float theta;
			float ro_dot;
			iss >> ro;
			iss >> theta;
			iss >> ro_dot;
			meas_package.raw_measurements_ << ro, theta, ro_dot;
			iss >> timestamp;
			meas_package.timestamp_ = timestamp;
			measurement_pack_list.push_back(meas_package);
		}
		//read gt
		float x_gt;
		float y_gt;
		float vx_gt;
		float vy_gt;
		iss >> x_gt;
		iss >> y_gt;
		iss >> vx_gt;
		iss >> vy_gt;
		gt_package.gt_values_ = VectorXd(4);
		gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
		gt_pack_list.push_back(gt_package);

	}

	//Create a UKF instance
	UKF ukf;

	//Call the UKF-based fusion
	size_t N = measurement_pack_list.size();
	for (size_t k = 0; k < N; ++k) {//start filtering from the second frame (the speed is unknown in the first frame)
	  ukf.ProcessMeasurement(measurement_pack_list[k]);

		//output the estimation
		out_file_ << ukf.x_(0) << "\t";
		out_file_ << ukf.x_(1) << "\t";
		out_file_ << ukf.x_(2) << "\t";
		out_file_ << ukf.x_(3) << "\t";

		//output the measurements
		if (measurement_pack_list[k].sensor_type_
				== MeasurementPackage::LASER) {
			//output the estimation
			out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
			out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
		} else if (measurement_pack_list[k].sensor_type_
				== MeasurementPackage::RADAR) {
			//output the estimation in the cartesian coordinates
			float ro = measurement_pack_list[k].raw_measurements_(0);
			float phi = measurement_pack_list[k].raw_measurements_(1);
			out_file_ << ro * cos(phi) << "\t";
			out_file_ << ro * sin(phi) << "\t";
		}

		//output the gt packages
		out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
		out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
		out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
		out_file_ << gt_pack_list[k].gt_values_(3) << "\n";

	}

	//close files
	if (out_file_.is_open()) {
		out_file_.close();
	}

	if (in_file.is_open()) {
		in_file.close();
	}
	return 0;
}

