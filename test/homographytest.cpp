/*
 * usacmethod.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#include "homographytest.h"
#include <map>
#include <algorithm>

vector<unsigned int> HomographyTest::match(const vector<cv::Point2f> &ref, const vector<cv::Point2f> &img, const vector< vector<int> > &inputCorrespondences)
{

	//Uniform sampling
	int maxCorresNum = inputCorrespondences.size();
	pointData.reserve(maxCorresNum);

	for (size_t i = 0 ; i < maxCorresNum ; i++)
	{
		pointData.push_back(ref[ inputCorrespondences[i][0] ].x);
		pointData.push_back(ref[ inputCorrespondences[i][0] ].y);
		pointData.push_back(1.0);
		pointData.push_back(img[ inputCorrespondences[i][1] ].x);
		pointData.push_back(img[ inputCorrespondences[i][1] ].y);
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

	homog = new usac::HomogEstimator;
	homog->initParamsUSAC(cfg);

	cfg.common.numDataPoints = maxCorresNum;


	//Create data index for prosac
	//The more reliable a correspondences, the higher the weight
	//We need to put correspondence with high reliability at the beginning of dataIndex array.

	vector<int> quality(maxCorresNum);
	for (size_t i = 0 ; i < maxCorresNum ; i++)
	{
		quality[i] = inputCorrespondences[i][2];
	}

	dataIndex.resize(maxCorresNum);
	std::vector<int>::iterator p;
	for (int i = 0 ; i < maxCorresNum ; i++)
	{
		p = std::max_element(quality.begin(),quality.end());
		dataIndex[i] = std::distance(quality.begin(), p);
		*p = -1;
	}

	cfg.prosac.sortedPointIndices = &dataIndex[0];


	/**
	 * Start to solve!
	 */
	homog->initDataUSAC(cfg);
	homog->initProblem(cfg, &pointData[0]);
	if (!homog->solve())
	{
		return vector<unsigned int>();     //no homo found
	}

	if (homog->usac_results_.best_inlier_count_ <= 8)
	{
		cv::Mat_<float> emptymat(3,3);
		emptymat.release();
		return vector<unsigned int>();     //no homo found
	}

	return homog->usac_results_.inlier_flags_;
}


/**
 * Read point set coordinate from specified file
 */
template <typename PointType>
bool HomographyTest::readCoordinateFromFile(string filename, vector<PointType> &pointarray)
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
bool HomographyTest::readResultFromFile(string filename, map<int,int> &resultmap)
{
	ifstream file(filename.c_str());
	if (!file.is_open())         //Check if the corresponding file exists
	{
		cout<<"Fail to open file "<< filename << endl;
		return false;
	}
	// Deal with all the coordinates

	int recordnum, refID, imgID;

	file >> recordnum;

	resultmap.clear();
	for(int i=0;i<recordnum;i++)
	{
		file >> refID >> imgID;
		resultmap[refID] = imgID;
	}

	file.close();

	return true;
}


bool HomographyTest::verify(string refFilename, string imgFilename, string resFilename)
{
	vector<cv::Point2f> refPts, imgPts;
	map<int,int> groundTruth;
	readCoordinateFromFile("/home/liming/workspace/usac/test/data/random50.ref", refPts);
	readCoordinateFromFile("/home/liming/workspace/usac/test/data/random50.img", imgPts);
	readResultFromFile("/home/liming/workspace/usac/test/data/random50.res", groundTruth);

	/**
	 * Matrix recording the existing correspondences
	 */
	vector< vector<bool> > isFilled(refPts.size());
	for (size_t i = 0 ; i < isFilled.size() ; i++)
	{
		isFilled[i].assign(imgPts.size() , false);
	}

	/**
	 * Input data (refPt.x refPt.y imgPt.x imgPt.y quality)
	 */
	vector< vector<int> > inputData (0);
	for (map<int,int>::iterator it = groundTruth.begin() ; it != groundTruth.end() ; ++it)
	{
		isFilled[it->first][it->second] = true;
		vector<int> onedata = {it->first, it->second, 2};
		inputData.push_back(onedata);
	}


	/**
	 * Then add noise from random data
	 */
	for (int i = 0 ; i < 2*groundTruth.size() ; i++)
	{
		int refID = int (rand() % refPts.size());
		int imgID = int (rand() % imgPts.size());
		if (!isFilled[refID][imgID])
		{
			isFilled[refID][imgID] = true;
			vector<int> onedata = {refID, imgID, 1};
			inputData.push_back(onedata);
		}
	}

	std::random_shuffle ( inputData.begin(), inputData.end() );
	vector<unsigned int> result = match(refPts, imgPts, inputData);

	map<int,int> resultMap;
	for (int i = 0 ; i < result.size() ; i++)
	{
		//Check the result here
		if (result[i])
		{
			resultMap[ inputData[i][0] ] = inputData[i][1];
		}
	}

	return (resultMap == groundTruth);
}

