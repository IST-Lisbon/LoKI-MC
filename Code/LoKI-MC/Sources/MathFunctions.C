#include "LoKI-MC/Headers/MathFunctions.h"
#include "LoKI-MC/Headers/Message.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <vector>
#include <iostream>
#include <cmath>
#include <random>

// used in the parallelization of the random number generation
// taken from https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
#include <thread>
#if defined (_MSC_VER)  // Visual studio
    #define thread_local __declspec( thread )
#else
    #define thread_local __thread
#endif

Eigen::ArrayXd MathFunctions::vectorToArray(std::vector<double> &vec){
	int vecSize = vec.size();
	Eigen::ArrayXd array(vec.size());
	for (int i = 0; i < vecSize; ++i){
		array[i] = vec[i];
	}
	return array;
}

double MathFunctions::unitUniformRand(bool includeZero, bool includeOne){
	// 'unitUniformRand' generates a random number according with a uniform distribution between 0 and 1
	// this function is thread-safe. See https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
	static thread_local std::random_device rd;
	static thread_local std::mt19937_64 generator(rd());
	static thread_local std::uniform_real_distribution<double> unitUniformDist(0.0, 1.0);
	double unitRandomNumber = unitUniformDist(generator);
	// loop to avoid zeros and ones, depending on the input parameters of the function
	while ((!includeZero && unitRandomNumber == 0) || (!includeOne && unitRandomNumber == 1)){
		unitRandomNumber = unitUniformDist(generator);
	}
	return unitRandomNumber;
}

double MathFunctions::unitNormalRand(){
	// 'unitNormalRand' generates a random number according with a normal distribution with mean value 0 and width 1
	// this function is thread-safe. See https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
	static thread_local std::random_device rd;
	static thread_local std::mt19937_64 generator(rd());
	static thread_local std::normal_distribution<double> unitNormalDist(0.0, 1.0);
	return unitNormalDist(generator);
}

Eigen::ArrayXd MathFunctions::histogramCount(Eigen::ArrayXd &arrayToCount, Eigen::ArrayXd &referenceGrid){
	// 'histogramCount' returns a cell histogram, given the nodes of the reference grid and the array elements to be counted
  	//
  	//  Node values -> |     |     |     |     |     |     |     |
	//  Cell values ->    o     o     o     o     o     o     o

	int numberOfCells = referenceGrid.size()-1;
	int arrayToCountSize = arrayToCount.size();
	double firstNode = referenceGrid[0];
	double step = referenceGrid[1]-firstNode;
	Eigen::ArrayXd arrayWithCounts = Eigen::ArrayXd::Zero(numberOfCells);
	// loop along all the elements that will be counted
	for (int i = 0; i < arrayToCountSize; ++i){
		int cellIndex = (arrayToCount[i]-firstNode)/step;
		if (cellIndex < numberOfCells){
			++arrayWithCounts[cellIndex];
		}
	}

	return arrayWithCounts;
}

Eigen::ArrayXd MathFunctions::histogramWeightedCount(Eigen::ArrayXd &arrayToCount, Eigen::ArrayXd &arrayToCountWeights, Eigen::ArrayXd &referenceGrid){
	// 'histogramWeightedCount' returns a cell histogram, given the nodes of the reference grid and the array elements to be counted 
  	//  Additionaly, each element can be counted with a different weight, using arrayToCountWeights
  	//
  	//  Node values -> |     |     |     |     |     |     |     |
	//  Cell values ->    o     o     o     o     o     o     o

	int numberOfCells = referenceGrid.size()-1;
	int arrayToCountSize = arrayToCount.size();
	double firstNode = referenceGrid[0];
	double step = referenceGrid[1]-firstNode;
	Eigen::ArrayXd arrayWithCounts = Eigen::ArrayXd::Zero(numberOfCells);
	// loop along all the elements that will be counted
	for (int i = 0; i < arrayToCountSize; ++i){
		int cellIndex = (arrayToCount[i]-firstNode)/step;
		if (cellIndex >= 0 && cellIndex < numberOfCells){
			arrayWithCounts[cellIndex] += arrayToCountWeights[i];
		}
	}

	return arrayWithCounts;
}

Eigen::ArrayXXd MathFunctions::histogram2DCount(Eigen::ArrayXd arrayToCountX, Eigen::ArrayXd arrayToCountY, Eigen::ArrayXd referenceGridX, Eigen::ArrayXd referenceGridY){
	// 'histogram2DCount' returns a 2D cell histogram, given the nodes of the 2 reference grids and the array elements to be counted

	int numberOfCellsX = referenceGridX.size()-1, numberOfCellsY = referenceGridY.size()-1;
	int arraysToCountSize = arrayToCountX.size();
	if (arraysToCountSize != arrayToCountY.size()){
		Message::error("The 2 arrays to be counted in 'histogram2DCount' should have equal dimensions!");
	}
	double firstNodeX = referenceGridX[0], firstNodeY = referenceGridY[0];
	double stepX = referenceGridX[1]-firstNodeX, stepY = referenceGridY[1]-firstNodeY;

	Eigen::ArrayXXd matrixWithCounts = Eigen::ArrayXXd::Zero(numberOfCellsX, numberOfCellsY);
	// loop along all the elements that will be counted
	for (int i = 0; i < arraysToCountSize; ++i){
		int cellIndexX = (arrayToCountX[i]-firstNodeX)/stepX;
		int cellIndexY = (arrayToCountY[i]-firstNodeY)/stepY;
		if (cellIndexX < numberOfCellsX && cellIndexX >= 0 && cellIndexY < numberOfCellsY && cellIndexY >= 0){
			++matrixWithCounts(cellIndexX,cellIndexY);
		}
	}
	return matrixWithCounts;
}

Eigen::Array3d MathFunctions::eulerTransformation(double chi, double eta, double theta, double phi){
	// according with Yousfi 1994	
	Eigen::Array3d versorInLab;
	double sinChi = std::sin(chi), cosChi = std::cos(chi), sinEta = std::sin(eta), cosEta = std::cos(eta), 
	       sinTheta = std::sin(theta), cosTheta = std::cos(theta), sinPhi = std::sin(phi), cosPhi = std::cos(phi);
	
	versorInLab[0] = -sinChi*sinEta*sinPhi + sinChi*cosEta*cosTheta*cosPhi + cosChi*sinTheta*cosPhi;
	versorInLab[1] = sinChi*sinEta*cosPhi + sinChi*cosEta*cosTheta*sinPhi + cosChi*sinTheta*sinPhi;
	versorInLab[2] = -sinChi*cosEta*sinTheta + cosChi*cosTheta;

	return versorInLab;
}

void MathFunctions::cart2sph(Eigen::Array3d &cartesianArray, double &norm, double &polarAngle, double &azimutalAngle){
	// 'cart2sph' calculates the spherical coordinates of an array
	// polarAngle: [0,pi]; azimutal angle [0,2pi]
	double x = cartesianArray[0], y = cartesianArray[1], z = cartesianArray[2];
	norm = std::sqrt(x*x + y*y + z*z);
	if (norm == 0){
		Message::error("Trying to use the function 'cart2sph' for a vector with null norm!");
	}
	azimutalAngle = std::atan2(y,x);
	// convert from [-pi,pi] to [0,2pi]
	if (azimutalAngle < 0){
		azimutalAngle += 2.0*M_PI;
	} 
	polarAngle = std::acos(z/norm);
}

Eigen::ArrayXd MathFunctions::standardDeviationColwise(Eigen::ArrayXXd matrix){
	int colSize = matrix.cols();
	Eigen::ArrayXd stds(colSize);
	for (int i = 0; i < colSize; ++i){
		stds[i] = MathFunctions::standardDeviation(matrix.col(i));
	}
	return stds;
}

double MathFunctions::standardDeviation(Eigen::ArrayXd array){
	double meanValue = array.mean();
	if (array.size() > 1){
		return (array - meanValue).matrix().norm() / std::sqrt(array.size()-1.0);
	}
	else{
		return 0;
	}
}

Eigen::ArrayXd MathFunctions::statisticalErrorColwise(Eigen::ArrayXXd matrix, double nBins){
	// 'statisticalErrorColwise' performs 'statisticalError' for each column

	int colSize = matrix.cols();
	Eigen::ArrayXd stds(colSize);
	for (int i = 0; i < colSize; ++i){
		stds[i] = statisticalError(matrix.col(i), nBins);
	}
	return stds;
}

double MathFunctions::statisticalError(Eigen::ArrayXd array, double nBins){
	// 'statisticalError' calculates the statistical error of the average of an array of 'nElements', dividing the 'nElements' in 'nBins' groups. 
	// This is done since each point may not be statistically independent, so we group them to get a realistic error estimation
	// Typically, nBins = 50

	int nElements = array.size();
	if (nElements == 0){
		return 0;
	}
	else if (nElements < nBins){
		Message::error("Error in 'statisticalError'. nBins is larger than nElements.");
	}

	int nPoints = nElements/nBins;
	double meanValue = array.mean();
	double sum = 0;
	for (int i = 0; i < nBins; ++i){
		sum += std::pow(array.segment(i*nPoints,nPoints).mean()-meanValue, 2);
	}
	return std::sqrt(sum)/nBins;
}


