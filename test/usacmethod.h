/*
 * usacmethod.h
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#ifndef USACLIB_USACMETHOD_H_
#define USACLIB_USACMETHOD_H_

#include "../basiclib/basictype.h"
#include "ConfigParamsHomog.h"
#include "HomogEstimator.h"
#include "USAC.h"

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
	PtArray refpts;
	PtArray imgpts;
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
	void setRefPoints(PtArray const *ptsource);

	/**
	 * Set reference point sets
	 */
	bool setImgPoints(PtArray const *ptsource);

	/**
	 * Matching with depth-first-search
	 */
	cv::Mat_<float> setImgPointsAndMatch(vector<RDM_Point> const *ptsource, const vector< int >* correspdist, vector<pair<int, pair<RDM_Point,RDM_Point> > > *corresp);

	/**
	 * Matching
	 */
	cv::Mat_<float> Match(vector<pair<int, pair<RDM_Point,RDM_Point> > > *corresp);

};

}

#endif /* USACLIB_USACMETHOD_H_ */
