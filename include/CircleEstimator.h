#ifndef CircleEstimator_H
#define CircleEstimator_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <Eigen/Dense>
#include "MathFunctions.h"
#include "CircleFittingFunctions.h"
#include "USAC.h"

namespace usac
{

class ConfigParamsCircle: public ConfigParams
{
	//To be edited if new parameters should be added
public:
	ConfigParamsCircle() {};
	virtual ~ConfigParamsCircle() {};
};

class CircleEstimator: public USAC<CircleEstimator>
{
	public:
		inline bool		 initProblem(const ConfigParamsCircle& cfg, double fixPoint[3], double* pointData);

		// ------------------------------------------------------------------------
		// problem specific functions: implement these
		inline void		 cleanupProblem();
		inline unsigned int generateMinimalSampleModels();
		inline bool		 generateRefinedModel(std::vector<unsigned int>& sample, const unsigned int numPoints, 
										  bool weighted = false, double* weights = NULL);
		inline bool		 validateSample();
		inline bool		 validateModel(unsigned int modelIndex);
		inline bool		 evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested);
		inline void		 testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel);
		inline unsigned int upgradeDegenerateModel();
		inline void		 findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
									 unsigned int numInliers, double* weights);
		inline void		 storeModel(unsigned int modelIndex, unsigned int numInliers);

	public:
		// ------------------------------------------------------------------------
		// storage for the final result
		std::vector<double> final_model_params_;


	private:
		double*      input_points_denorm_;
		double 	 fix_point_[3], fix_point_denorm_[3];     //The circle should pass this fix point
		// ------------------------------------------------------------------------
		// temporary storage
		double* input_points_;							// stores normalized data points
//		double* data_matrix_;							// linearized input data
		double  m_T_[9], m_Tinv_[9];					// normalization matrices
		std::vector<double*> models_;				    // stores vector of models, (x0, y0,r)
		std::vector<double*> models_denorm_;			// stores vector of (denormalized) models
};


// ============================================================================================
// initProblem: initializes problem specific data and parameters
// call this function once per run on new data
// pointData = 3*n double, (in homogeneous cooridnates) [x1, y1, 1, x2, y2, 1, x3, y3, 1...]
// ============================================================================================
bool CircleEstimator::initProblem(const ConfigParamsCircle& cfg, double fixPoint[3], double* pointData)
{
	// set up problem specific data here (e.g., data matrices, parameter vectors, etc.)
	std::memcpy(fix_point_denorm_, fixPoint, sizeof(fix_point_));

	// copy pointer to input data
	input_points_denorm_ = pointData;
	input_points_        = new double[6*cfg.common.numDataPoints];
	if (input_points_denorm_ == NULL)
	{
		std::cerr << "Input point data not properly initialized" << std::endl;
		return false;
	}
	if (input_points_ == NULL)
	{
		std::cerr << "Could not allocate storage for normalized data points" << std::endl;
		return false;
	}

	// normalize input data
	// following this, input_points_ has the normalized points and input_points_denorm_ has
	// the original input points
	// output: input_points_
	MathTools::normalizePoints(input_points_denorm_, input_points_, cfg.common.numDataPoints, m_T_);
	MathTools::vmul( fix_point_, m_T_, fix_point_denorm_, 3);


	for (unsigned int i = 0; i < 9; ++i)
	{
		m_Tinv_[i] = m_T_[i];
	}
	MathTools::minv(m_Tinv_, 3);

	// allocate storage for models
	size_t model_parameter_number = 3;
	final_model_params_.clear(); final_model_params_.resize(model_parameter_number);
	models_.clear(); models_.resize(usac_max_solns_per_sample_);
	models_denorm_.clear(); models_denorm_.resize(usac_max_solns_per_sample_);
	for (unsigned int i = 0; i < usac_max_solns_per_sample_; ++i)
	{
		models_[i] = new double[model_parameter_number];
		models_denorm_[i] = new double[model_parameter_number];
	}

	// precompute the data matrix
//	data_matrix_ = new double[18*usac_num_data_points_];	// 2 equations per correspondence
//	HTools::computeDataMatrix(data_matrix_, usac_num_data_points_, input_points_);

	return true;
}


// ============================================================================================
// cleanupProblem: release any temporary problem specific data storage 
// call this function once at the end of each run on new data
// ============================================================================================
void CircleEstimator::cleanupProblem()
{
	if (input_points_) { delete[] input_points_; input_points_ = NULL; }
//	if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
	for (size_t i = 0; i < models_.size(); ++i)
	{
		if (models_[i]) { delete[] models_[i]; }
	}
	models_.clear();
	for (size_t i = 0; i < models_denorm_.size(); ++i)
	{
		if (models_denorm_[i]) { delete[] models_denorm_[i]; }
	}
	models_denorm_.clear();
}


// ============================================================================================
// generateMinimalSampleModels: generates minimum sample model from the data points whose
// indices are currently stored in "min_sample_[]".
// in this case, only one model per minimum sample
// ============================================================================================
unsigned int CircleEstimator::generateMinimalSampleModels()
{
	// this function must be implemented
	// Step 3 in USAC::solve()

	std::vector<Eigen::Vector2f> sourcePoints(usac_min_sample_size_);
	for (unsigned int i = 0; i < usac_min_sample_size_; ++i)
	{
		sourcePoints[i](0) = *(input_points_ + 3*min_sample_[i] + 0);
		sourcePoints[i](1) = *(input_points_ + 3*min_sample_[i] + 1);
	}

	sourcePoints.push_back( Eigen::Vector2f(float(fix_point_[0]), float(fix_point_[1])  ) );

	Eigen::Vector2f center;
	float radius;
	CircleTools::circleFit(sourcePoints, center, radius);
	*(models_[0] + 2) = radius;
	*(models_[0] + 0) = center(0);   *(models_[0] + 1) = center(1);

	double centerPos[3] = {center(0), center(1), 1};
	double centerPos_denorm[3];

	MathTools::vmul( centerPos_denorm, m_Tinv_, centerPos, 3);
	*(models_denorm_[0] + 2) = radius*m_Tinv_[0];
	*(models_denorm_[0] + 0) = centerPos_denorm[0];   *(models_denorm_[0] + 1) = centerPos_denorm[1];


	return 1;  // return the number of minimal sample models
}


// ============================================================================================
// generateRefinedModel: compute model using non-minimal set of samples
// default operation is to use a weight of 1 for every data point
// ============================================================================================
bool CircleEstimator::generateRefinedModel(std::vector<unsigned int>& sample,
										  unsigned int numPoints,
										  bool weighted,
										  double* weights)
{
	// implement this function if you want to use local optimization
	std::cout<<"This is called only with locallyOptimizeSolution."<<std::endl;
	return true;
}


// ============================================================================================
// validateSample: check if minimal sample is valid
// ============================================================================================
bool CircleEstimator::validateSample()
{
   // implement this function if you want to pre-verify the minimal sample, otherwise
   // simply return true

   //Step 2 in USAC::solve(), for circle fitting, any sample is OK.
   return true;	
}


// ============================================================================================
// validateModel: check if model computed from minimal sample is valid
// ============================================================================================
bool CircleEstimator::validateModel(const unsigned int modelIndex)
{
	// Step 4 in USAC::solve()
   // implement this function if you want to pre-verify a model (controled by usac_prevalidate_model_), otherwise
   // simply return true
	return true;
}


// ============================================================================================
// evaluateModel: test model against all/subset of the data points
// ============================================================================================
bool CircleEstimator::evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested)
{
	// test model against all data points, or a randomized subset (in the case of early
	// termination with, for e.g., SPRT)
	double* model = models_denorm_[modelIndex];
	std::vector<double>::iterator current_err_array = err_ptr_[0];
	double err;
	bool good_flag = true;
	double lambdaj, lambdaj_1 = 1.0;
	*numInliers = 0;
	*numPointsTested = 0;
	unsigned int pt_index;

	//Added
	double* pt;

	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		// get index of point to be verified
		if (eval_pool_index_ > usac_num_data_points_-1)
		{
			eval_pool_index_ = 0;
		}
		pt_index = evaluation_pool_[eval_pool_index_];
		++eval_pool_index_;

		// compute point-model error for the problem of interest
		//
		// --- implement this
		// pt_index is the index of the point to be verified
		// Get the distance from the point to the circle as the error
		//
		// Step 5 in USAC::solve()

		pt = input_points_denorm_ + 3*pt_index;   //get the pointer to the point
		double dx = pt[0] - model[0];
		double dy = pt[1] - model[1];
		err = fabs( sqrt(dx*dx + dy*dy) - model[2] );
		//end of computer error

		*(current_err_array+pt_index) = err;

		if (err < usac_inlier_threshold_)
		{
			++(*numInliers);
		}

		if (usac_verif_method_ == USACConfig::VERIF_SPRT)
		{
			if (err < usac_inlier_threshold_)
			{			
				lambdaj = lambdaj_1 * (sprt_delta_/sprt_epsilon_);
			}
			else
			{
				lambdaj = lambdaj_1 * ( (1 - sprt_delta_)/(1 - sprt_epsilon_) );
			}

			if (lambdaj > decision_threshold_sprt_)
			{
				good_flag = false;
				*numPointsTested = i+1;
				return good_flag;
			}
			else
			{
				lambdaj_1 = lambdaj;
			}
		}
	}
	*numPointsTested = usac_num_data_points_;

	return good_flag;
}

// ============================================================================================
// testSolutionDegeneracy: check if model is degenerate
// ============================================================================================
void CircleEstimator::testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel)
{
	// implement this if you want to detect degenerate models, otherwise return false
	*degenerateModel = false;
	*upgradeModel = false;
}

// ============================================================================================
// upgradeDegenerateModel: try to upgrade degenerate model to non-degenerate by sampling from
// the set of outliers to the degenerate model
// ============================================================================================
unsigned int CircleEstimator::upgradeDegenerateModel()
{
	// implement this if you want to upgrade degenerate models to non-degenerate, otherwise
	// return 0;
	return 0;
}


// ============================================================================================
// findWeights: given model and points, compute weights to be used in local optimization
// ============================================================================================
void CircleEstimator::findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers,
								 unsigned int numInliers, double* weights)
{
	// implement this if you want weighted local refinement, otherwise return an array of ones
	for (unsigned int i = 0; i < numInliers; ++i)
	{
		weights[i] = 1.0;
	}
}


// ============================================================================================
// storeModel: stores current best model
// this function is called  (by USAC) every time a new best model is found
// ============================================================================================
void CircleEstimator::storeModel(const unsigned int modelIndex, unsigned int numInliers)
{
	// save the current model as the best solution so far
	for (unsigned int i = 0; i < 3; ++i)
	{
		final_model_params_[i] = *(models_denorm_[modelIndex]+i);
	}
}

}

#endif

