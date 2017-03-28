/*
 * usacmethod.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#include "circleFittest.h"
#include <map>
#include <algorithm>

vector<unsigned int> CircleFitTest::fit(const cv::Point2f &fixPt, const vector<cv::Point2f> &pts, const vector< int > &q)
{

	//Uniform sampling
	int maxPoints = pts.size();
	pointData.reserve(maxPoints);

	for (size_t i = 0 ; i < maxPoints ; i++)
	{
		pointData.push_back(pts[ i ].x);
		pointData.push_back(pts[ i ].y);
		pointData.push_back(1.0);
	}

	//4 unit of noise in reference dataset
	//	jitter = 4;

	//In comment, we keep the default parameters
//	cfg.common.minSampleSize = 4;
//	cfg.common.maxSolutionsPerSample = 1;
//	cfg.common.maxHypotheses = 100000;
	cfg.common.randomSamplingMethod = USACConfig::SAMP_PROSAC;
	cfg.common.inlierThreshold = 2*sqrt(2)*jitter;   //2*sqrt(2)*sigma, for symetric error

	//Optimizer
	cfg.common.prevalidateSample = true;
	cfg.common.prevalidateModel = true;
	cfg.common.testDegeneracy = true;
	cfg.common.localOptMethod = USACConfig::LO_LOSAC;

	circleFittor = new usac::CircleEstimator;
	circleFittor->initParamsUSAC(cfg);

	cfg.common.numDataPoints = maxPoints;


	//Create data index for prosac
	//The more reliable a correspondences, the higher the weight
	//We need to put correspondence with high reliability at the beginning of dataIndex array.

	vector<int> quality(maxPoints);
	for (size_t i = 0 ; i < maxPoints ; i++)
	{
		quality[i] = q[i];
	}

	dataIndex.resize(maxPoints);
	std::vector<int>::iterator p;
	for (int i = 0 ; i < maxPoints ; i++)
	{
		p = std::max_element(quality.begin(),quality.end());
		dataIndex[i] = std::distance(quality.begin(), p);
		*p = -1;
	}

	cfg.prosac.sortedPointIndices = &dataIndex[0];


	/**
	 * Start to solve!
	 */
	double fix_point_xy1[3];
	fix_point_xy1[0] = fixPt.x; fix_point_xy1[1] = fixPt.y; fix_point_xy1[2] = 1;
	circleFittor->initDataUSAC(cfg);
	circleFittor->initProblem(cfg, fix_point_xy1, &pointData[0]);
	if (!circleFittor->solve())
	{
		return vector<unsigned int>();     //no homo found
	}

	if (circleFittor->usac_results_.best_inlier_count_ <= 8)
	{
		cv::Mat_<float> emptymat(3,3);
		emptymat.release();
		return vector<unsigned int>();     //no homo found
	}

	return circleFittor->usac_results_.inlier_flags_;
}


/**
 * Read point set coordinate from specified file
 */
template <typename PointType>
bool CircleFitTest::readPointsFromFile(string filename, vector<PointType> &pointarray)
{
	ifstream file(filename.c_str());
	if (!file.is_open())         //Check if the corresponding file exists
	{cv::Mat_<float> Match(vector<pair<int, pair<cv::Point2f,cv::Point2f> > > *corresp);
		cout<<"Fail to open file "<< filename << endl;
		return false;
	}
	// Deal with all the coordinates

	int totalpoints;
	PointType tmppoint;
	pointarray.clear();

	file >> totalpoints;

	for(int i=0;i<totalpoints;i++)
	{
		file >> tmppoint.x >> tmppoint.y;
		pointarray.push_back(tmppoint);
	}

	file.close();

	return true;
}

/**
 * Read expected result from specified file, the result is recorded in a map
 * resultmap[refID] = imgID;
 */
bool CircleFitTest::readResultFromFile(string filename, vector<int> &resultmap)
{
	ifstream file(filename.c_str());
	if (!file.is_open())         //Check if the corresponding file exists
	{
		cout<<"Fail to open file "<< filename << endl;
		return false;
	}
	// Deal with all the coordinates

	int recordnum, ptID;

	file >> recordnum;

	resultmap.clear();
	for(int i=0;i<recordnum;i++)
	{
		file >> ptID;
		resultmap.push_back(ptID);
	}

	file.close();

	return true;
}


bool CircleFitTest::verify(string refFilename, string imgFilename, string resFilename)
{
	vector<cv::Point2f> pts;
	vector<int> groundTruth;
	readPointsFromFile("/home/liming/workspace/usac/test/data/random50.ref", pts);
	readResultFromFile("/home/liming/workspace/usac/test/data/random50.res", groundTruth);

}

