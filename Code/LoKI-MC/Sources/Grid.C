#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>

Grid::Grid(){

	std::string eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");
	if (eedfType == "prescribedEedf"){
		cellNumber = FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.cellNumber");
		double maxEnergy = FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.maxEnergy");	
		step = maxEnergy/cellNumber;
		node = Eigen::ArrayXd::LinSpaced(cellNumber+1, 0, maxEnergy);
		cell = Eigen::ArrayXd::LinSpaced(cellNumber, step/2.0, maxEnergy-step/2.0);

		if (FieldInfo::getField("electronKinetics.numerics.energyGrid.smartGrid")){
			isSmart = true;
			minEedfDecay = FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay");
			maxEedfDecay = FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay");
			updateFactor = FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.updateFactor");
		}
	}
	else{
		if (FieldInfo::getField("electronKinetics.numericsMC.nEnergyCells")){
			cellNumber = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nEnergyCells");
		}
		else{
			cellNumber = 500;
		}
		double maxEnergy = 1; // this is irrelevant at this stage, since the maxEnergy is defined in the Monte Carlo simulation
		step = maxEnergy/cellNumber;
		node = Eigen::ArrayXd::LinSpaced(cellNumber+1, 0, maxEnergy);
		cell = Eigen::ArrayXd::LinSpaced(cellNumber, step/2.0, maxEnergy-step/2.0);
	}
}

void Grid::updateMaxValue(double maxValue){
	// resize the grid with a new maximum value
	step = maxValue/cellNumber;
	node = Eigen::ArrayXd::LinSpaced(cellNumber+1, 0, maxValue);
	cell = Eigen::ArrayXd::LinSpaced(cellNumber, step/2.0, maxValue-step/2.0);

	// change whatever is necessary due to the update of the max energy value
	updatedMaxEnergy1Signal();
	updatedMaxEnergy2Signal();
}

Eigen::ArrayXd Grid::adjustToGrid(std::vector<std::vector<double>> originalArray, double threshold, std::string mode){
	Eigen::ArrayXd x;
	if (mode == "nodes"){
		x = node;
	}
	else if (mode == "cells"){
		x = cell;
	}
	else{
		Message::error(std::string("Unrecognized mode '") + mode + " while adjusting data to Grid\n");
	}

	int xSize = x.size(), originalArraySize = originalArray[0].size(), i;
	Eigen::ArrayXd adjustedArray = Eigen::ArrayXd::Zero(xSize);
	double x1, x2, y1, y2;
	if (threshold < x[xSize-1]){
		i = 1;
		for (int j = ceil(threshold/step); j<xSize; ++j){
			x1 = originalArray[0][i-1];
			x2 = originalArray[0][i];
			while (x2<x[j] && i<originalArraySize){
				++i;
				x1 = x2;
				x2 = originalArray[0][i];
			}
			if (x2 >= x[j]){
				y1 = originalArray[1][i-1];
				y2 = originalArray[1][i];
				adjustedArray[j] = y1 + (x[j]-x1)*(y2-y1)/(x2-x1);
				if (adjustedArray[j]<0){
					adjustedArray[j] = 0;
				}
			}
			else{
				adjustedArray[j] = 0;
			}
		}
	}

	return adjustedArray;
}

double Grid::integrate(Eigen::ArrayXd integrand){
	// integrate evaluates the integral of the array "integrand" with respect to the Grid 'grid' using the trapezoidal rule
	double result = 0;
	int integrandSize = integrand.size();
	if (integrandSize == cell.size()){
		for (int i = 0; i < integrandSize; ++i){
			result += integrand[i] * step;
		}
	}
	else if (integrandSize == node.size()){
		for (int i = 0; i < integrandSize-1; ++i){
			result += step*0.5*(integrand[i]+integrand[i+1]);		
		}
	}
	else{
		Message::error("Error! Integrand dimensions do not match with grid dimensions.\n");
	}
	return result;
}

double Grid::maxwellianRateCoeff(Eigen::ArrayXd crossSection, double temperatureInEV){
    // maxwellianRateCoeff evaluates the rate coefficient of a certain collision with cross section 'crossSection' whith a maxwellian eedf
    // characterised by temperature 'temperatureInEV'. The 'crossSection' array must be adjusted to the 'grid', and the 'temperatureInEV' must
    // be in electron volts.
    // result = sqrt(2e/m_e)*integral(sigma(u)*u*f(u)*du)

	Eigen::ArrayXd f;
	double result = 0;

    if (crossSection.size() == cell.size()){
    	f = Eigen::exp(cell/(-temperatureInEV));
    	f /= integrate( f*Eigen::sqrt(cell) );
    	result = integrate(crossSection*cell*f);
    }
    else if (crossSection.size() == node.size()){
    	f = Eigen::exp(node/(-temperatureInEV));
    	f /= integrate( f*Eigen::sqrt(node) );
    	result = integrate(crossSection*node*f);
    }
    else{
    	Message::error("Integrand dimensions do not match with grid dimensions.");
    }   

    result *= std::sqrt(2.0*Constant::electronCharge/Constant::electronMass);

    return result;
}