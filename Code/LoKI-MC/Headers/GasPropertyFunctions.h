#ifndef __GasPropertyFunctions__
#define __GasPropertyFunctions__
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Constant.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

//****************************************************//
//   property functions related with heat capacity    //
//****************************************************//

template <class GasType>
double nitrogenHeatCapacity(GasType* gas, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	// heat capacity for nitrogen (SI units)
	double heatCapacity = 29.1 + 2494.2/(553.4*std::sqrt(M_PI/2.0))*std::exp(-2.0*std::pow((workCond->gasTemperature-1047.4)/553.4, 2));
	// change energy units to eV
	heatCapacity /= Constant::electronCharge;
	// store value on gas object properties
	gas->heatCapacity = heatCapacity;
	return heatCapacity;	
}

template <class GasType>
double oxygenHeatCapacity(GasType* gas, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	// heat capacity for oxygen (SI units)
	double heatCapacity = 28.8 + 6456.2/(788.3*std::sqrt(M_PI/2.0))*std::exp(-2.0*std::pow((workCond->gasTemperature-1006.9)/788.3, 2));
	// change energy units to eV
	heatCapacity /= Constant::electronCharge;
	// store value on gas object properties
	gas->heatCapacity = heatCapacity;
	return heatCapacity;	
}

//****************************************************//
//property functions related with thermal conductivity//
//****************************************************//

template <class GasType>
double nitrogenThermalConductivity(GasType* gas, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	// thermal conductivity for nitrogen (SI units)
	double thermalConductivity = (1.717+0.084*workCond->gasTemperature-1.948e-5*std::pow(workCond->gasTemperature,2))*1e-3;
	// change energy units to eV
	thermalConductivity /= Constant::electronCharge;
	// store value on gas object properties
	gas->thermalConductivity = thermalConductivity;
	return thermalConductivity;
}

template <class GasType>
double oxygenThermalConductivity(GasType* gas, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	// thermal conductivity for oxygen (SI units)
	double thermalConductivity = (1.056+0.087*workCond->gasTemperature-8.912e-6*std::pow(workCond->gasTemperature,2))*1e-3;
	// change energy units to eV
	thermalConductivity /= Constant::electronCharge;
	// store value on gas object properties
	gas->thermalConductivity = thermalConductivity;
	return thermalConductivity;
}

template <class GasType>
void constantValue(GasType* gas, std::string property, std::vector<double> &argumentArray){
	double value = argumentArray[0];
	if (property == "mass"){
		gas->mass = value;
	}
	else if (property == "harmonicFrequency"){
		gas->harmonicFrequency = value;
	}
	else if (property == "anharmonicFrequency"){
		gas->anharmonicFrequency = value;
	}
	else if (property == "rotationalConstant"){
		gas->rotationalConstant = value;
	}
	else if (property == "lennardJonesDistance"){
		gas->lennardJonesDistance = value;
	}
	else if (property == "lennardJonesDepth"){
		gas->lennardJonesDepth = value;
	}
	else if (property == "electricDipolarMoment"){
		gas->electricDipolarMoment = value;
	}
	else if (property == "electricQuadrupoleMoment"){
		gas->electricQuadrupoleMoment = value;
	}
	else if (property == "polarizability"){
		gas->polarizability = value;
	}
	else if (property == "fraction"){
		gas->fraction = value;
	}
	else if (property == "heatCapacity"){
		gas->heatCapacity = value;
	}
	else if (property == "thermalConductivity"){
		gas->thermalConductivity = value;
	}
	else if (property == "OPBParameter"){
		gas->OPBParameter = value;
	}
	else{
		Message::error(std::string("Error. The gas property '") + property + "' is not recognized by the code.\n Please check the input file.");
	}
}

template <class GasType>
void evaluateGasPropertyFunction(std::string propertyString, GasType* gas, std::string property, std::vector<double> argumentArray, std::vector<std::string> argumentStrArray, WorkingConditions* workCond){
	// 'evaluateGasPropertyFunction' evaluates the property function 'propertyFuncStr' for a given gas.
	// The corresponding property function pointer is saved in the gas, to be used later on if necessary
	//
	// add here the new property functions
	
	std::string propertyFuncStr = Parse::tokenizeCharacters(propertyString,(char*) "@")[0];

	// ----- property functions related with heat capacity ----- //

	if (propertyFuncStr == "nitrogenHeatCapacity"){
		// check if the parameters are appropriate to this function
		if (property != "heatCapacity"){
			Message::error(std::string("Trying to use nitrogenHeatCapacity function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating nitrogenHeatCapacity function. \nCheck input file");
		}
		// calculate the values of the property
		nitrogenHeatCapacity(gas, argumentArray, argumentStrArray, workCond);
		// assign this function to the gas (to be used later on, if necessary)
		gas->heatCapacityFunc = nitrogenHeatCapacity;
		gas->heatCapacityParams = argumentArray;
		gas->heatCapacityStrParams = argumentStrArray;
	}

	else if (propertyFuncStr == "oxygenHeatCapacity"){
		// check if the parameters are appropriate to this function
		if (property != "heatCapacity"){
			Message::error(std::string("Trying to use oxygenHeatCapacity function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating oxygenHeatCapacity function. \nCheck input file");
		}
		// calculate the values of the property
		oxygenHeatCapacity(gas, argumentArray, argumentStrArray, workCond);
		// assign this function to the gas (to be used later on, if necessary)
		gas->heatCapacityFunc = oxygenHeatCapacity;
		gas->heatCapacityParams = argumentArray;
		gas->heatCapacityStrParams = argumentStrArray;
	}

	// ----- property functions related with thermal conductivity ----- //

	else if (propertyFuncStr == "nitrogenThermalConductivity"){
		// check if the parameters are appropriate to this function
		if (property != "thermalConductivity"){
			Message::error(std::string("Trying to use nitrogenThermalConductivity function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating nitrogenThermalConductivity function. \nCheck input file");
		}
		// calculate the values of the property
		nitrogenThermalConductivity(gas, argumentArray, argumentStrArray, workCond);
		// assign this function to the gas (to be used later on, if necessary)
		gas->thermalConductivityFunc = nitrogenThermalConductivity;
		gas->thermalConductivityParams = argumentArray;
		gas->thermalConductivityStrParams = argumentStrArray;
	}

	else if (propertyFuncStr == "oxygenThermalConductivity"){
		// check if the parameters are appropriate to this function
		if (property != "thermalConductivity"){
			Message::error(std::string("Trying to use oxygenThermalConductivity function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating oxygenThermalConductivity function. \nCheck input file");
		}
		// calculate the values of the property
		oxygenThermalConductivity(gas, argumentArray, argumentStrArray, workCond);
		// assign this function to the gas (to be used later on, if necessary)
		gas->thermalConductivityFunc = oxygenThermalConductivity;
		gas->thermalConductivityParams = argumentArray;
		gas->thermalConductivityStrParams = argumentStrArray;
	}

	else if (propertyString.find("@") != std::string::npos){
		Message::error(std::string("Error! Trying to use the property function '") + propertyFuncStr + "' which is not defined in the code.\n" +
			           "Add it at 'GasPropertyFunctions.h', and do not forget to put it also at the function 'evaluateGasPropertyFunction'.");
	}

	else{
		argumentArray.push_back(Parse::str2value(propertyString));
		constantValue(gas, property, argumentArray);
	}
}

#endif