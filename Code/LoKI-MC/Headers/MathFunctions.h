#ifndef __MathFunctions__
#define __MathFunctions__

#include "External/eigen-3.4.0/Eigen/Dense"
#include <vector>

// eigen array of booleans
namespace Eigen{
	typedef Array<bool,Dynamic,1> ArrayXb;
};

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
	Eigen::Array3d unitNormalRand3();
	Eigen::ArrayXd histogramCount(Eigen::ArrayXd &arrayToCount, Eigen::ArrayXd &referenceGrid);
	Eigen::ArrayXd histogramWeightedCount(Eigen::ArrayXd &arrayToCount, Eigen::ArrayXd &arrayToCountWeights, Eigen::ArrayXd &referenceGrid);
	Eigen::ArrayXXd histogram2DCount(Eigen::ArrayXd arrayToCountX, Eigen::ArrayXd arrayToCountY, Eigen::ArrayXd referenceGridX, Eigen::ArrayXd referenceGridY);
	Eigen::Array3d eulerTransformation(double sinChi, double cosChi, double sinEta, double cosEta, double sinTheta, double cosTheta, double sinPhi, double cosPhi);
	void cart2sph(Eigen::Array3d &cartesianArray, double &norm, double &sinPolarAngle, double &cosPolarAngle, double &sinAzimutalAngle, double &cosAzimutalAngle);

	Eigen::ArrayXd standardDeviationColwise(Eigen::ArrayXXd matrix);
	double standardDeviation(Eigen::ArrayXd array);
	Eigen::ArrayXd statisticalErrorColwise(Eigen::ArrayXXd matrix, double nPoints);
	double statisticalError(Eigen::ArrayXd array, double nPoints);
};

#endif