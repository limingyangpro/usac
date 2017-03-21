/*
 * usacmethod.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: liming
 */

#include "usacmethod.h"
#include <map>

namespace usac
{
void USACmethod::setJitter()
{
	int num,i;
	float refarea;                   //minimum rectangle englobe all the referece points
	cv::Point2f minp,maxp;             //minimum rectangle englobe all the referece points
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

void USACmethod::setRefPoints(vector<cv::Point2f> const *ptsource)
{
	refpts = *ptsource;
	setJitter();
}

/**
 * Set reference point sets
 */
bool USACmethod::setImgPoints(vector<cv::Point2f> const *ptsource)
{
	//Simple copy make core dump, why ?
	imgpts.resize(ptsource->size());
	for (int i = 0 ; i < ptsource->size() ; i++)
	{
		imgpts[i] = (*ptsource)[i];
	}
	return true;
}

cv::Mat_<float> USACmethod::Match(vector<pair<int, pair<cv::Point2f,cv::Point2f> > > *corresp)
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

cv::Mat_<float> USACmethod::setImgPointsAndMatch(vector<cv::Point2f> const *ptsource, const vector< int >* correspdist, vector<pair<int, pair<cv::Point2f,cv::Point2f> > > *corresp)
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


vector<unsigned int> USACmethod::Match(const vector<cv::Point2f> &ref, const vector<cv::Point2f> &img, const vector< vector<int> > &inputCorrespondences)
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

	jitter = 4;
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


}

/**
 * Read point set coordinate from specified file
 */
template <typename PointType>
bool readCoordinateFromFile(string filename, vector<PointType> &pointarray)
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
bool readResultFromFile(string filename, map<int,int> &resultmap)
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


int main()
{
	vector<cv::Point2f> refPts, imgPts;
	map<int,int> groundTruth;

	readCoordinateFromFile("/home/liming/workspace/usac/test/data/grafitti1and2.ref", refPts);
	readCoordinateFromFile("/home/liming/workspace/usac/test/data/grafitti1and2.img", imgPts);
	readResultFromFile("/home/liming/workspace/usac/test/data/grafitti1and2.res", groundTruth);

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
		vector<int> onedata = {it->first, it->second, 1};
		inputData.push_back(onedata);
	}


	/**
	 * Then add noise from random data
	 */
//	for (int i = 0 ; i < 2*groundTruth.size() ; i++)
//	{
//		int refID = int (rand() % refPts.size());
//		int imgID = int (rand() % imgPts.size());
//		if (!isFilled[refID][imgID])
//		{
//			isFilled[refID][imgID] = true;
//			vector<int> onedata = {refID, imgID, 1};
//			inputData.push_back(onedata);
//		}
//	}

	usac::USACmethod usacalgo(0.05, 1);

	vector<unsigned int> result = usacalgo.Match(refPts, imgPts, inputData);

	for (int i = 0 ; i < result.size() ; i++)
	{
		cout<<result[i]<<endl;
	}

	return 0;

}
