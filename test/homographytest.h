/*
 * usacmethod.h
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#ifndef HOMOGRAPHYTEST_H_
#define HOMOGRAPHYTEST_H_

#include "ConfigParamsHomog.h"
#include "HomogEstimator.h"
#include "USAC.h"
#include "opencv2/opencv.hpp"
using namespace std;

struct Corresp
{
	int refID, imgID;
	int quality;
	Corresp(int ref, int img, int quality_) : refID(ref), imgID(img), quality(quality_) {};
};

class HomographyTest
{
private:
	float jitter;

	std::vector<double> pointData;
	std::vector<unsigned int> dataIndex;
	usac::ConfigParamsHomog cfg;
	usac::HomogEstimator* homog;

public:
	HomographyTest(float jitter_) :
		jitter(jitter_){};

	/**
	 * New method for match
	 */
	vector<unsigned int> match(const vector<cv::Point2f> &ref, const vector<cv::Point2f> &img, const vector< vector<int> > &inputCorrespondences);

	template <typename PointType>
	bool readCoordinateFromFile(string filename, vector<PointType> &pointarray);
	bool readResultFromFile(string filename, map<int,int> &resultmap);

	bool verify(string refFilename, string imgFilename, string resFilename);

};

#endif /* HOMOGRAPHYTEST_H_ */
