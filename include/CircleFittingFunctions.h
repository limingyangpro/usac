#ifndef CTOOLS_H
#define CTOOLS_H

#include "MathFunctions.h"
#include <Eigen/Dense>

namespace CircleTools
{
//	void computeDataMatrix(double* data_matrix, unsigned int num_points, double* points);
	void circleFit(const std::vector<Eigen::Vector2f>& sourcePoints, Eigen::Vector2f &center, float &radius);
//	void crossprod(double *out, const double *a, const double *b, unsigned int st);
}
#endif
