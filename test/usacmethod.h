/*
 * usacmethod.h
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#ifndef USACLIB_USACMETHOD_H_
#define USACLIB_USACMETHOD_H_

#include "ConfigParamsHomog.h"
#include "HomogEstimator.h"
#include "USAC.h"
#include "opencv2/opencv.hpp"
using namespace std;

namespace usac
{

struct Corresp
{
	int refID, imgID;
	int quality;
	Corresp(int ref, int img, int quality_) : refID(ref), imgID(img), quality(quality_) {};
};

class USACmethod
{
private:
	vector<cv::Point2f> refpts;
	vector<cv::Point2f> imgpts;
	float eta, jitter;
	int distThreshold;
	const vector< int >* correspondancedist;  /**< refID*imgnum + imgID: quality of correspondence, from "color distance", the smaller the better */

	std::vector<double> pointData;
	std::vector<unsigned int> dataIndex;
	ConfigParamsHomog cfg;
	HomogEstimator* homog;

private:
	void setJitter();
	bool setCorrespondenceDist(const vector< int >* correspdist);


public:
	USACmethod(float eta_, int distThreshold_) :
		eta(eta_),
		distThreshold(distThreshold_),
		jitter(-1),
		correspondancedist(NULL),
		homog(NULL) {};

	/**
	 * Set reference point sets and renew the hashtables
	 */
	void setRefPoints(vector<cv::Point2f> const *ptsource);

	/**
	 * Set reference point sets
	 */
	bool setImgPoints(vector<cv::Point2f> const *ptsource);

	/**
	 * Matching with depth-first-search
	 */
	cv::Mat_<float> setImgPointsAndMatch(vector<cv::Point2f> const *ptsource, const vector< int >* correspdist, vector<pair<int, pair<cv::Point2f,cv::Point2f> > > *corresp);

	/**
	 * Matching
	 */
	cv::Mat_<float> Match(vector<pair<int, pair<cv::Point2f,cv::Point2f> > > *corresp);

	/**
	 * New method for match
	 */
	vector<unsigned int> Match(const vector<cv::Point2f> &ref, const vector<cv::Point2f> &img, const vector< vector<int> > &inputCorrespondences);
};

}

#endif /* USACLIB_USACMETHOD_H_ */
