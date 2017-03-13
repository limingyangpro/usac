/*
 * usacmethod.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#include "usacmethod.h"

namespace usac
{
void USACmethod::setJitter()
{
	int num,i;
	float refarea;                   //minimum rectangle englobe all the referece points
	RDM_Point minp,maxp;             //minimum rectangle englobe all the referece points
	num = refpts.size();

	//Find bounding area
	minp = refpts[0];
	maxp = refpts[0];
	for(i = 1 ; i < num ; i++)
	{
		//Note this can not work in 3D
		if (refpts[i].x < minp.x) minp.x = refpts[i].x;
		if (refpts[i].x > maxp.x) maxp.x = refpts[i].x;
		if (refpts[i].y < minp.y) minp.y = refpts[i].y;
		if (refpts[i].y > maxp.y) maxp.y = refpts[i].y;
	}
	refarea = (maxp.x - minp.x)*(maxp.y - minp.y);

	jitter = eta*sqrt(refarea/num);
}

bool USACmethod::setCorrespondenceDist(const vector< int >* correspdist)
{
	correspondancedist = correspdist;

	if (correspondancedist != NULL)
	{
		//Find max dist
		int maxDist = *std::max_element(correspondancedist->begin(),correspondancedist->end());

		//create valid correspondence vector
		int imgNum = imgpts.size();
		int refNum = refpts.size();
		int maxCorresNum = imgNum*refNum*0.7;

		std::vector<int>LCvotes;
		LCvotes.clear();
		LCvotes.reserve(maxCorresNum);
		pointData.reserve(maxCorresNum*6);

		for (int r = 0 ; r < refNum ; r++)
		{
			for (int i = 0 ; i < imgNum ; i++)
				if ((*correspondancedist)[r*imgNum + i] < distThreshold)
				{
					pointData.push_back(refpts[r].x);
					pointData.push_back(refpts[r].y);
					pointData.push_back(1.0);
					pointData.push_back(imgpts[i].x);
					pointData.push_back(imgpts[i].y);
					pointData.push_back(1.0);
					LCvotes.push_back(maxDist - (*correspondancedist)[r*imgNum + i]);
				}
		}

		//Sort the quality into dataindex
		dataIndex.resize(LCvotes.size());
		std::vector<int>::iterator p;
		for (int i = 0 ; i < LCvotes.size() ; i++)
		{
			p = std::max_element(LCvotes.begin(),LCvotes.end());
			dataIndex[i] = std::distance(LCvotes.begin(), p);
			*p = -1;
		}

//		for (int i = 0 ; i < LCvotes.size() ; i++)
//		{
//			cout<<LCvotes[i]<<" order: "<<dataIndex[i]<<endl;
//		}

		cfg.common.randomSamplingMethod = USACConfig::SAMP_PROSAC;
	//	cfg.common.minSampleSize = 4;
	//	cfg.common.maxSolutionsPerSample = 1;
	//	cfg.common.maxHypotheses = 100000;
	//  cfg.common.randomSamplingMethod = USACConfig::SAMP_UNIFORM;
		cfg.common.inlierThreshold = 2*sqrt(2)*jitter;   //2*sqrt(2)*sigma, for symetric error

		//param for prosac
		double		  beta;
		double        nonRandConf;
		//Optimizer
		cfg.common.prevalidateSample = true;
		cfg.common.prevalidateModel = true;
		cfg.common.testDegeneracy = true;

		cfg.common.localOptMethod = USACConfig::LO_LOSAC;

		homog = new HomogEstimator;
		homog->initParamsUSAC(cfg);

		cfg.common.numDataPoints = LCvotes.size();

		cfg.prosac.sortedPointIndices = &dataIndex[0];
	}
	else
	{
		//create valid correspondence vector
		int imgNum = imgpts.size();
		int refNum = refpts.size();
		int maxCorresNum = imgNum*refNum;
		pointData.reserve(maxCorresNum*6);

		for (int r = 0 ; r < refNum ; r++)
		{
			for (int i = 0 ; i < imgNum ; i++)
			{
				pointData.push_back(refpts[r].x);
				pointData.push_back(refpts[r].y);
				pointData.push_back(1.0);
				pointData.push_back(imgpts[i].x);
				pointData.push_back(imgpts[i].y);
				pointData.push_back(1.0);
			}
		}

	//	cfg.common.minSampleSize = 4;
	//	cfg.common.maxSolutionsPerSample = 1;
	//	cfg.common.maxHypotheses = 100000;
		cfg.common.randomSamplingMethod = USACConfig::SAMP_UNIFORM;
		cfg.common.inlierThreshold = 2*sqrt(2)*jitter;   //2*sqrt(2)*sigma, for symetric error
		cfg.common.prevalidateSample = true;
		cfg.common.prevalidateModel = true;
		cfg.common.testDegeneracy = true;

		cfg.common.localOptMethod = USACConfig::LO_LOSAC;

		homog = new HomogEstimator;
		homog->initParamsUSAC(cfg);

		cfg.common.numDataPoints = maxCorresNum;
	}
}

void USACmethod::setRefPoints(PtArray const *ptsource)
{
	refpts = *ptsource;
	setJitter();
}

/**
 * Set reference point sets
 */
bool USACmethod::setImgPoints(PtArray const *ptsource)
{
	//Simple copy make core dump, why ?
	imgpts.resize(ptsource->size());
	for (int i = 0 ; i < ptsource->size() ; i++)
	{
		imgpts[i] = (*ptsource)[i];
	}
	return true;
}

cv::Mat_<float> USACmethod::Match(vector<pair<int, pair<RDM_Point,RDM_Point> > > *corresp)
{
	cv::Mat_<float> homo = cv::Mat::eye(3,3,CV_32F);
	// set up the homography estimation problem
	homog->initDataUSAC(cfg);
	homog->initProblem(cfg, &pointData[0]);
	if (!homog->solve())
	{
		cv::Mat_<float> emptymat(3,3);
		emptymat.release();
		return emptymat;     //no homo found
	}

	if (homog->usac_results_.best_inlier_count_ <= 8)
	{
		cv::Mat_<float> emptymat(3,3);
		emptymat.release();
		return emptymat;     //no homo found
	}

	//output homography
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			homo(i,j) = (float)homog->final_model_params_[3*i+j]/homog->final_model_params_[8];
		}
	}

	return homo;
}

cv::Mat_<float> USACmethod::setImgPointsAndMatch(vector<RDM_Point> const *ptsource, const vector< int >* correspdist, vector<pair<int, pair<RDM_Point,RDM_Point> > > *corresp)
{
	if (setImgPoints(ptsource))
	{
		//Set usac lib parameters
		setCorrespondenceDist(correspdist);

		return Match(corresp);
	}
	cv::Mat_<float> emptymat(3,3);
	emptymat.release();
	return emptymat;
}


}
