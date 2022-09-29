#ifndef __MathFunctions__
#define __MathFunctions__

#include "External/eigen-3.4.0/Eigen/Dense"
#include <vector>

namespace MathFunctions{
	Eigen::ArrayXd vectorToArray(std::vector<double> &vec);
	template <class ObjectType>
	std::vector<ObjectType> append(std::vector<ObjectType> v1, std::vector<ObjectType> v2){
		for (auto value: v2){
			v1.push_back(value);
		}
		return v1;
	}


	double unitUniformRand(bool includeZero, bool includeOne);
	double unitNormalRand();
	Eigen::ArrayXd histogramCount(Eigen::ArrayXd &arrayToCount, Eigen::ArrayXd &referenceGrid);
	Eigen::ArrayXd histogramWeightedCount(Eigen::ArrayXd &arrayToCount, Eigen::ArrayXd &arrayToCountWeights, Eigen::ArrayXd &referenceGrid);
	Eigen::ArrayXXd histogram2DCount(Eigen::ArrayXd arrayToCountX, Eigen::ArrayXd arrayToCountY, Eigen::ArrayXd referenceGridX, Eigen::ArrayXd referenceGridY);
	Eigen::Array3d eulerTransformation(double chi, double eta, double theta, double phi);
	void cart2sph(Eigen::Array3d &cartesianArray, double &norm, double &polarAngle, double &azimutalAngle);

	Eigen::ArrayXd standardDeviationColwise(Eigen::ArrayXXd matrix);
	double standardDeviation(Eigen::ArrayXd array);
	Eigen::ArrayXd statisticalErrorColwise(Eigen::ArrayXXd matrix, double nPoints);
	double statisticalError(Eigen::ArrayXd array, double nPoints);
};

#endif