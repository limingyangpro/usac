/*
 * usacmethod.h
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#ifndef USACLIB_USACMETHOD_H_
#define USACLIB_USACMETHOD_H_

#include "CircleEstimator.h"
#include "USAC.h"
#include "opencv2/opencv.hpp"
using namespace std;

class CircleFitTest
{
private:
	float jitter;

	std::vector<double> pointData;
	std::vector<unsigned int> dataIndex;
	usac::ConfigParamsCircle cfg;
	usac::CircleEstimator* circleFittor;

public:
	CircleFitTest(float jitter_) :
		jitter(jitter_){};

	/**
	 * New method for match
	 */
	vector<unsigned int> fit(const cv::Point2f &fixPt, const vector<cv::Point2f> &pts, const vector< int > &quatlity);

	template <typename PointType>
	bool readPointsFromFile(string filename, vector<PointType> &pointarray);
	bool readResultFromFile(string filename, vector<int> &resultmap);

	bool verify(string refFilename, string imgFilename, string resFilename);

};

#endif /* USACLIB_USACMETHOD_H_ */
