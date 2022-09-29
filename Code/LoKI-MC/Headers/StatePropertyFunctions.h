#ifndef __StatePropertyFunctions__
#define __StatePropertyFunctions__
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/WorkingConditions.h"
#include <string>
#include <vector>
#include <cmath>

//****************************************************//
// property functions related with statistical weights//
//****************************************************//

template <class StateType>
void rotationalDegeneracy(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double J;
	for (auto& state: stateArray){
		J = std::stod(state->rotLevel);
		state->statisticalWeight = 2.0*J + 1.0;
	}
}

template <class StateType>
void rotationalDegeneracy_H2(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){	
	double J;
	for (auto& state: stateArray){
		J = std::stod(state->rotLevel);
		state->statisticalWeight = (2.0-std::pow(-1,J))*(2.0*J+1.0);
	}
}

template <class StateType>
void rotationalDegeneracy_N2(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double J;
	for (auto& state: stateArray){
		J = std::stod(state->rotLevel);
		state->statisticalWeight = 3.0*(1.0+0.5*(1.0+std::pow(-1,J)))*(2.0*J+1.0);
	}
}

template <class StateType>
void rotationalDegeneracy_H2O(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
    // the rot. state is given as a string JKaKc
    // split the string into individual components
    for (auto& state: stateArray){
        double J = std::stod(state->rotLevel.substr(0,1));
        double Ka = std::stod(state->rotLevel.substr(1,1));
        double Kc = std::stod(state->rotLevel.substr(2,1));
        // statistical weight
        double s = 2.0*J + 1;
        // distinguish ortho and para states
        double tau = std::abs(Ka-Kc);
        if (std::fmod(tau,2) != 0){
        	s *= 3.0;
        }
        state->statisticalWeight = s;
    }
}
//****************************************************//
// property functions related with state energies	  //
//****************************************************//

template <class StateType>
void harmonicOscillatorEnergy(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double vibLevel;
	for (auto& state: stateArray){
		vibLevel = std::stod(state->vibLevel);
		state->energy = Constant::planckReducedInEV * state->gas->harmonicFrequency*(vibLevel+0.5);
	}
}

template <class StateType>
void morseOscillatorEnergy(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double vibLevel;
	for(auto& state: stateArray){
		vibLevel = std::stod(state->vibLevel);
		state->energy = Constant::planckReducedInEV * (state->gas->harmonicFrequency*(vibLevel+0.5)-state->gas->anharmonicFrequency*std::pow((vibLevel+0.5),2));
	}
}

template <class StateType>
void rigidRotorEnergy(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double rotLevel;
	for(auto& state: stateArray){
		rotLevel = std::stod(state->rotLevel);
		state->energy = state->gas->rotationalConstant*rotLevel*(rotLevel+1.0);
	}
}


//****************************************************//
// property functions related with diffusion/mobility //
//****************************************************//

template <class StateType>
void generalizedTemperatureDependentDiffCoeff(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double normalizationConstant = argumentArray[0];
	double temperature;
	if (argumentStrArray[1] == "gasTemperature"){
		temperature = workCond->gasTemperature;
	}
	else if (argumentStrArray[1] == "electronTemperature"){
		temperature = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temperature = argumentArray[1];
	}
	double power = argumentArray[2];
	for (auto& state: stateArray){
		state->reducedDiffCoeff = normalizationConstant*std::pow(temperature,power);
	}		
}

template <class StateType>
void generalizedTemperatureDependentMobCoeff(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double normalizationConstant = argumentArray[0];
	double temperature;
	if (argumentStrArray[1] == "gasTemperature"){
		temperature = workCond->gasTemperature;
	}
	else if (argumentStrArray[1] == "electronTemperature"){
		temperature = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temperature = argumentArray[1];
	}
	double power = argumentArray[2];
	for (auto& state: stateArray){
		state->reducedMobility = normalizationConstant*std::pow(temperature,power);
	}
}

//****************************************************//
// property functions related with state populations  //
//****************************************************//

template <class StateType>
void boltzmannPopulation(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double temperature;
	if (argumentStrArray[0] == "gasTemperature"){
		temperature = workCond->gasTemperature;
	}
	else if (argumentStrArray[0] == "electronTemperature"){
		temperature = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temperature = argumentArray[0];
	}
	double norm = 0;
	// calculate the ground state energy. This is important to avoid problems when the temperature is too low
	double groundStateEnergy = 1E20;
	for (auto& state: stateArray){
		groundStateEnergy = std::fmin(groundStateEnergy, state->energy);
	}
	for (auto& state: stateArray){
		state->population = state->statisticalWeight * std::exp(-(state->energy-groundStateEnergy)/(Constant::boltzmannInEV*temperature));
		norm += state->population;
	}
	for (auto& state: stateArray){
		state->population /= norm;
		state->populationIsUpdated = true;
	}	
}

template <class StateType>
void boltzmannPopulationVibrationalCutoff(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double temperature;
	if (argumentStrArray[0] == "gasTemperature"){
		temperature = workCond->gasTemperature;
	}
	else if (argumentStrArray[0] == "electronTemperature"){
		temperature = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temperature = argumentArray[0];
	}
	// obtain the cutoff vibrational level for the Boltzmann distribution (this level is included in the distribution)
	double vMax = argumentArray[1];
	double norm = 0;
	// calculate the ground state energy. This is important to avoid problems when the temperature is too low
	double groundStateEnergy = 1E20;
	for (auto& state: stateArray){
		groundStateEnergy = std::fmin(groundStateEnergy, state->energy);
	}
	for (auto& state: stateArray){
		if (std::stod(state->vibLevel) > vMax){
			state->population = 0;
		}
		else{
			state->population = state->statisticalWeight * std::exp(-(state->energy-groundStateEnergy)/(Constant::boltzmannInEV*temperature));
			norm += state->population;
		}
	}
	for (auto& state: stateArray){
		state->population /= norm;
		state->populationIsUpdated = true;
	}	
}

template <class StateType>
void boltzmannPopulationRotationalCutoff(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double temperature;
	if (argumentStrArray[0] == "gasTemperature"){
		temperature = workCond->gasTemperature;
	}
	else if (argumentStrArray[0] == "electronTemperature"){
		temperature = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temperature = argumentArray[0];
	}
	double jMax = argumentArray[1]; // cutoff rotational level for the Boltzmann distribution (this level is included in the distribution)
	double norm = 0;
	// calculate the ground state energy. This is important to avoid problems when the temperature is too low
	double groundStateEnergy = 1E20;
	for (auto& state: stateArray){
		groundStateEnergy = std::fmin(groundStateEnergy, state->energy);
	}
	for (auto& state: stateArray){
		if (std::stod(state->rotLevel) > jMax){
			state->population = 0;
		}
		else{
			state->population = state->statisticalWeight * std::exp(-(state->energy-groundStateEnergy)/(Constant::boltzmannInEV*temperature));
			norm += state->population;
		}
	}
	for (auto& state: stateArray){
		state->population /= norm;
		state->populationIsUpdated = true;
	}	
}

template <class StateType>
void treanorGordietsPopulation(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){
	double temp0;
	if (argumentStrArray[0] == "gasTemperature"){
		temp0 = workCond->gasTemperature;
	}
	else if (argumentStrArray[0] == "electronTemperature"){
		temp0 = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temp0 = argumentArray[0];
	}
	double temp1;
	if (argumentStrArray[1] == "gasTemperature"){
		temp1 = workCond->gasTemperature;
	}
	else if (argumentStrArray[1] == "electronTemperature"){
		temp1 = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temp1 = argumentArray[1];
	}
	double anharmonicFrequency = stateArray[0]->gas->anharmonicFrequency;

	double norm = 0, vibLevel, vLimit;
	double groundEnergy, firstEnergy;
	bool groundDefined = false, firstDefined = false;

	for (auto& state: stateArray){
		if (state->vibLevel == "0"){
			groundEnergy = state->energy;
			groundDefined = true;
		}
		else if (state->vibLevel == "1"){
			firstEnergy = state->energy;
			firstDefined = true;
		}
	}

	if ((!groundDefined) || (!firstDefined)){
		Message::error(std::string("Unable to find groundEnergy or firstEnergy to populate state") + stateArray[0]->name + " and its siblings "+
		"with function treanorGordietsPopulation. \nCheck input file.");
	}
	else{
		vLimit = std::floor( 0.5*(1.0+(firstEnergy-groundEnergy)*temp0/(Constant::planckReducedInEV*anharmonicFrequency*temp1)) );
	}

	// evaluate regular Treanor
	int stateArraySize = stateArray.size();
	int indexVLimit;
	StateType* state;
	for (int i = 0; i < stateArraySize; ++i){
		state = stateArray[i];
		vibLevel = std::stod(state->vibLevel);
		if (vibLevel == vLimit){
			indexVLimit = i;
		}
		state->population = state->statisticalWeight * std::exp( -1.0/Constant::boltzmannInEV * (vibLevel*(firstEnergy-groundEnergy)*
														   (1.0/temp1-1.0/temp0) + (state->energy-groundEnergy)/temp0));
		norm += state->population;
	}

	for (auto& state: stateArray){
		state->population /= norm;
	}

	// evaluate Treanor-Gordiets (modified Treanor above vLimit)
	std::vector<double> populationTreanorGordiets(stateArraySize,0);
	norm = 0;
	for (int i = 0; i<stateArraySize; ++i){
		state = stateArray[i];
		vibLevel = std::stod(state->vibLevel);
		if (vibLevel <= vLimit){
			populationTreanorGordiets[i] = state->population;
		}
		else{
			populationTreanorGordiets[i] = stateArray[indexVLimit]->population*vLimit/vibLevel;
		}
		norm += populationTreanorGordiets[i]; 
	}

	for (int i = 0; i<stateArraySize; ++i){
		stateArray[i]->population = populationTreanorGordiets[i]/norm;
		stateArray[i]->populationIsUpdated = true;
	}
}

template <class StateType>
void treanorPopulation(std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond){

	double temp0;
	if (argumentStrArray[0] == "gasTemperature"){
		temp0 = workCond->gasTemperature;
	}
	else if (argumentStrArray[0] == "electronTemperature"){
		temp0 = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temp0 = argumentArray[0];
	}
	double temp1;
	if (argumentStrArray[1] == "gasTemperature"){
		temp1 = workCond->gasTemperature;
	}
	else if (argumentStrArray[1] == "electronTemperature"){
		temp1 = workCond->electronTemperature*Constant::electronCharge/Constant::boltzmann;
	}
	else{
		temp1 = argumentArray[1];
	}
	double norm = 0, vibLevel;
	double groundEnergy, firstEnergy;
	bool groundDefined = false, firstDefined = false;

	for (auto state: stateArray){
		if (state->vibLevel == "0"){
			groundEnergy = state->energy;
			groundDefined = true;
		}
		else if (state->vibLevel == "1"){
			firstEnergy = state->energy;
			firstDefined = true;
		}
	}

	if ((!groundDefined) || (!firstDefined)){
		Message::error(std::string("Unable to find groundEnergy or firstEnergy to populate state") + stateArray[0]->name + " and its siblings "+
		"with function treanorPopulation. \nCheck input file.");
	}

	for (auto& state: stateArray){
		vibLevel = std::stod(state->vibLevel);
		state->population = state->statisticalWeight * std::exp( -1.0/Constant::boltzmannInEV * (vibLevel*(firstEnergy-groundEnergy)*
														   (1.0/temp1-1.0/temp0) + (state->energy-groundEnergy)/temp0));
		norm += state->population;
	}

	for (auto& state: stateArray){
		state->population /= norm;
		state->populationIsUpdated = true;
	}
}

template <class StateType>
void constantValue(std::vector<StateType*> &stateArray, std::string property, std::vector<double> &argumentArray){
	double value = argumentArray[0];
	if (property == "energy"){
		for (auto& state: stateArray){
			state->energy = value;
		}
	}
	else if (property == "statisticalWeight"){
		for (auto& state: stateArray){
			state->statisticalWeight = value;
		}
	}
	else if (property == "reducedDiffCoeff"){
		for (auto& state: stateArray){
			state->reducedDiffCoeff = value;
		}
	}
	else if (property == "reducedMobility"){
		for (auto& state: stateArray){
			state->reducedMobility = value;
		}
	}
	else if (property == "population"){
		for (auto& state: stateArray){
			state->population = value;
		}
	}
	else{
		Message::error(std::string("Error. The state property '") + property + "' is not recognized by the code.\n Please check the input file.");
	}
}

template <class StateType>
void evaluateStatePropertyFunction(std::string propertyString, std::vector<StateType*> &stateArray, std::string property, std::vector<double> argumentArray, std::vector<std::string> argumentStrArray, WorkingConditions* workCond){
	// 'evaluateStatePropertyFunction' evaluates the property function 'propertyFuncStr' for a given set of states.
	// The corresponding property function pointer is saved in the states, to be used later on if necessary
	//
	// add here the new property functions
	
	std::string propertyFuncStr = Parse::tokenizeCharacters(propertyString,(char*) "@")[0];

	// ----- property functions related with statistical weights ----- //

	if (propertyFuncStr == "rotationalDegeneracy"){
		// check if the parameters are appropriate to this function
		if (property != "statisticalWeight"){
			Message::error(std::string("Trying to use rotationalDegeneracy function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating rotationalDegeneracy function. \nCheck input file");
		}
		else if (stateArray[0]->type != "rot"){
			Message::error(std::string("Trying to assign rotationalDegeneracy to non rotational state ") + stateArray[0]->name + ". \nCheck input file");
		}
		// calculate the values of the property
		rotationalDegeneracy(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->statisticalWeightFunc = rotationalDegeneracy;
			state->statisticalWeightParams = argumentArray;
			state->statisticalWeightStrParams = argumentStrArray;
		}
	}

	else if (propertyFuncStr == "rotationalDegeneracy_H2"){
		// check if the parameters are appropriate to this function
		if (property != "statisticalWeight"){
			Message::error("Trying to use rotationalDegeneracy_H2 function to set up property " + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating rotationalDegeneracy_H2 function. \nCheck input file");
		}
		else if (stateArray[0]->type != "rot"){
			Message::error(std::string("Trying to assign rotationalDegeneracy_H2 to non rotational state ") + stateArray[0]->name + ". \nCheck input file");
		}
		// calculate the values of the property
		rotationalDegeneracy_H2(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->statisticalWeightFunc = rotationalDegeneracy_H2;
			state->statisticalWeightParams = argumentArray;
			state->statisticalWeightStrParams = argumentStrArray;		
		}		
	}

	else if (propertyFuncStr == "rotationalDegeneracy_N2"){
		// check if the parameters are appropriate to this function
		if (property != "statisticalWeight"){
			Message::error(std::string("Trying to use rotationalDegeneracy_N2 function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating rotationalDegeneracy_N2 function. \nCheck input file");
		}
		else if (stateArray[0]->type != "rot"){
			Message::error(std::string("Trying to assign rotationalDegeneracy_N2 to non rotational state ") + stateArray[0]->name + ". \nCheck input file");
		}
		// calculate the values of the property
		rotationalDegeneracy_N2(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->statisticalWeightFunc = rotationalDegeneracy_N2;
			state->statisticalWeightParams = argumentArray;
			state->statisticalWeightStrParams = argumentStrArray;		
		}	
	}

	else if (propertyFuncStr == "rotationalDegeneracy_H2O"){
		// check if the parameters are appropriate to this function
		if (property != "statisticalWeight"){
			Message::error(std::string("Trying to use rotationalDegeneracy_H2O function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating rotationalDegeneracy_H2O function. \nCheck input file");
		}
		else if (stateArray[0]->type != "rot"){
			Message::error(std::string("Trying to assign rotationalDegeneracy_H2O to non rotational state ") + stateArray[0]->name + ". \nCheck input file");
		}
		// calculate the values of the property
		rotationalDegeneracy_H2O(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->statisticalWeightFunc = rotationalDegeneracy_H2O;
			state->statisticalWeightParams = argumentArray;
			state->statisticalWeightStrParams = argumentStrArray;		
		}	
	}

	// ----- property functions related with state energies ----- //

	else if (propertyFuncStr == "harmonicOscillatorEnergy"){
		// check if the parameters are appropriate to this function
		if (property != "energy"){
			Message::error(std::string("Trying to use harmonicOscillatorEnergy function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating harmonicOscillatorEnergy function. \nCheck input file");
		}
		else if (stateArray[0]->type != "vib"){
			Message::error(std::string("Trying to assign harmonicOscillatorEnergy to non vibrational state ") + stateArray[0]->name + ".\nCheck input file.");
		}
		else if (stateArray[0]->gas->harmonicFrequency == Constant::NON_DEF){
			Message::error(std::string("Unable to find harmonicFrequency to evaluate the energy of the state ") + stateArray[0]->name + " with function harmonicOscillatorEnergy.\nCheck input file");
		}
		// calculate the values of the property
		harmonicOscillatorEnergy(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->energyFunc = harmonicOscillatorEnergy;
			state->energyParams = argumentArray;
			state->energyStrParams = argumentStrArray;
		}
	}

	else if (propertyFuncStr == "morseOscillatorEnergy"){
		// check if the parameters are appropriate to this function
		StateType* firstState = stateArray[0];
		if (property != "energy"){
			Message::error(std::string("Trying to use morseOscillatorEnergy function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating morseOscillatorEnergy function. \nCheck input file");
		}
		else if (firstState->type != "vib"){
			Message::error(std::string("Trying to assign morseOscillatorEnergy to non vibrational state ") + firstState->name + "\nCheck input file.\n");
		}
		else if (firstState->gas->harmonicFrequency == Constant::NON_DEF){
			Message::error(std::string("Unable to find harmonicFrequency to evaluate the energy of the state ") + firstState->name + " with function morseOscillatorEnergy.\nCheck input file");
		}
		else if (firstState->gas->anharmonicFrequency == Constant::NON_DEF){
			Message::error(std::string("Unable to find anharmonicFrequency to evaluate the energy of the state ") + firstState->name + " with function morseOscillatorEnergy.\nCheck input file");
		}
		// calculate the values of the property
		morseOscillatorEnergy(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->energyFunc = morseOscillatorEnergy;
			state->energyParams = argumentArray;
			state->energyStrParams = argumentStrArray;
		}		
	}

	else if (propertyFuncStr == "rigidRotorEnergy"){
		// check if the parameters are appropriate to this function
		StateType* firstState = stateArray[0];
		if (property != "energy"){
			Message::error(std::string("Trying to use rigidRotorEnergy function to set up property ") + property + ". \nCheck input file");
		}
		else if (!argumentArray.empty()){
			Message::error("Wrong number of arguments when evaluating rigidRotorEnergy function. \nCheck input file");
		}
		else if (firstState->type != "rot"){
			Message::error(std::string("Trying to assign rigidRotorEnergy to non rotational state ") + firstState->name + "\nCheck input file.");
		}
		else if (firstState->gas->rotationalConstant == Constant::NON_DEF){
			Message::error(std::string("Unable to find rotationalConstant to evaluate the energy of the state ") + firstState->name + " with function rigidRotorEnergy.\nCheck input file");
		}
		// calculate the values of the property
		rigidRotorEnergy(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->energyFunc = rigidRotorEnergy;
			state->energyParams = argumentArray;
			state->energyStrParams = argumentStrArray;
		}	
	}

	// ----- property functions related with diffusion/mobility	----- //

	else if (propertyFuncStr == "generalizedTemperatureDependentCoeff"){
		// check if the parameters are appropriate to this function
		if ((property != "reducedDiffCoeff") && (property != "reducedMobility")){
			Message::error(std::string("Trying to use generalizedTemperatureDependentCoeff to set up property ") + property + ".\nCheck input file.");
		}
		else if (argumentArray.size() != 3){
			Message::error("Wrong number of arguments when evaluating generalizedTemperatureDependentCoeff function. \nCheck input file");	
		}
		else if (property == "reducedDiffCoeff"){
			// calculate the values of the property
			generalizedTemperatureDependentDiffCoeff(stateArray, argumentArray, argumentStrArray, workCond);
			// assign this function to each state (to be used later on, if necessary)
			for (auto& state: stateArray){
				state->reducedDiffCoeffFunc = generalizedTemperatureDependentDiffCoeff;
				state->reducedDiffCoeffParams = argumentArray;
				state->reducedDiffCoeffStrParams = argumentStrArray;
			}
		}
		else if (property == "reducedMobility"){
			// calculate the values of the property
			generalizedTemperatureDependentMobCoeff(stateArray, argumentArray, argumentStrArray, workCond);
			// assign this function to each state (to be used later on, if necessary)
			for (auto& state: stateArray){
				state->reducedMobilityFunc = generalizedTemperatureDependentMobCoeff;
				state->reducedMobilityParams = argumentArray;
				state->reducedMobilityStrParams = argumentStrArray;
			}
		}	
	}

	// ----- property functions related with state populations ----- //

	else if (propertyFuncStr == "boltzmannPopulation"){
		// check if the parameters are appropriate to this function
		if (property != "population"){
			Message::error(std::string("Trying to use boltzmannPopulation to set up property ") + property + ".\nCheck input file.");
		}
		else if (argumentArray.size() != 1){
			Message::error("Wrong number of arguments when evaluating boltzmannPopulation function. \nCheck input file");	
		}
		// check if the energies and statistical weights are defined
		for (auto& state: stateArray){
			if (state->energy == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " energy for the evaluation of boltzmannPopulation function. \nCheck input file");
			}
			else if (state->statisticalWeight == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " statistical weight for the evaluation of boltzmannPopulation function. \nCheck input file\n");
			}
		}
		// calculate the values of the property
		boltzmannPopulation(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->populationFunc = boltzmannPopulation;
			state->populationParams = argumentArray;
			state->populationStrParams = argumentStrArray;
		}		
	}

	else if (propertyFuncStr == "boltzmannPopulationVibrationalCutoff"){
		// check if the parameters are appropriate to this function
		if (property != "population"){
			Message::error(std::string("Trying to use boltzmannPopulationVibrationalCutoff to set up property ") + property + ".\nCheck input file.");
		}
		else if (argumentArray.size() != 2){
			Message::error("Wrong number of arguments when evaluating boltzmannPopulationVibrationalCutoff function. \nCheck input file");	
		}
		else if (stateArray[0]->type != "vib"){
			Message::error(std::string("Trying to assign boltzmannPopulationVibrationalCutoff to non vibrational state") + stateArray[0]->name + ". \nCheck input file");
		}
		double vMax = argumentArray[1];
		if (std::fmod(vMax,1) != 0 || vMax < 0){
			Message::error(std::string("Error found when evaluating population of state ") + stateArray[0]->name + ".\nMaximum vibrational level ''" + std::to_string(vMax) + "should be a non-negative integer.\nCheck input file");
		}
		// check if the energies and statistical weights are defined
		for (auto& state: stateArray){
			if (state->energy == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " energy for the evaluation of boltzmannPopulationVibrationalCutoff function. \nCheck input file");
			}
			else if (state->statisticalWeight == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " statistical weight for the evaluation of boltzmannPopulationVibrationalCutoff function. \nCheck input file\n");
			}
		}
		// calculate the values of the property
		boltzmannPopulationVibrationalCutoff(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->populationFunc = boltzmannPopulationVibrationalCutoff;
			state->populationParams = argumentArray;
			state->populationStrParams = argumentStrArray;
		}
	}

	else if (propertyFuncStr == "boltzmannPopulationRotationalCutoff"){
		// check if the parameters are appropriate to this function
		if (property != "population"){
			Message::error(std::string("Trying to use boltzmannPopulationRotationalCutoff to set up property ") + property + ".\nCheck input file.");
		}
		else if (argumentArray.size() != 2){
			Message::error("Wrong number of arguments when evaluating boltzmannPopulationRotationalCutoff function. \nCheck input file");	
		}
		else if (stateArray[0]->type != "rot"){
			Message::error(std::string("Trying to assign boltzmannPopulationRotationalCutoff to non rotational state") + stateArray[0]->name + ". \nCheck input file");
		}
		double jMax = argumentArray[1];
		if (std::fmod(jMax,1) != 0 || jMax < 0){
			Message::error(std::string("Error found when evaluating population of state ") + stateArray[0]->name + ".\nMaximum rotational level ''" + std::to_string(jMax) + "should be a non-negative integer.\nCheck input file");
		}
		// check if the energies and statistical weights are defined
		for (auto& state: stateArray){
			if (state->energy == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " energy for the evaluation of boltzmannPopulationRotationalCutoff function. \nCheck input file");
			}
			else if (state->statisticalWeight == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " statistical weight for the evaluation of boltzmannPopulationRotationalCutoff function. \nCheck input file\n");
			}
		}
		// calculate the values of the property
		boltzmannPopulationRotationalCutoff(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->populationFunc = boltzmannPopulationRotationalCutoff;
			state->populationParams = argumentArray;
			state->populationStrParams = argumentStrArray;
		}		
	}

	else if (propertyFuncStr == "treanorGordietsPopulation"){
		// check if the parameters are appropriate to this function
		if (property != "population"){
			Message::error(std::string("Trying to use treanorGordietsPopulation function to set up property ") + property + ". \nCheck input file");
		}
		else if (argumentArray.size() != 2){
			Message::error("Wrong number of arguments when evaluating treanorGordietsPopulation function. \nCheck input file");
		}
		else if (stateArray[0]->gas->anharmonicFrequency == Constant::NON_DEF){
			Message::error(std::string("Unable to find ") + stateArray[0]->gas->name + " anharmonicFrequency to populate" + stateArray[0]->name
			+ " and its siblings with function treanorGordietsPopulation. \nCheck input file\n");
		}
		else if (stateArray[0]->type != "vib"){
			Message::error(std::string("Trying to assign treanorGordietsPopulation to non vibrational state") + stateArray[0]->name + ". \nCheck input file");
		}
		// check if the energies and statistical weights are defined
		for (auto& state: stateArray){
			if (state->energy == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " energy for the evaluation of treanorGordietsPopulation function. \nCheck input file");
			}
			else if (state->statisticalWeight == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " statistical weight for the evaluation of treanorGordietsPopulation function. \nCheck input file");
			}
		}
		// calculate the values of the property
		treanorGordietsPopulation(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->populationFunc = treanorGordietsPopulation;
			state->populationParams = argumentArray;
			state->populationStrParams = argumentStrArray;
		}			
	}

	else if (propertyFuncStr == "treanorPopulation"){
		// check if the parameters are appropriate to this function
		if (property != "population"){
			Message::error(std::string("Trying to use treanorPopulation function to set up property ") + property + ". \nCheck input file");
		}
		else if (argumentArray.size() != 2){
			Message::error("Wrong number of arguments when evaluating treanorPopulation function. \nCheck input file");
		}
		else if (stateArray[0]->type != "vib"){
			Message::error(std::string("Trying to assign treanorPopulation to non vibrational state") + stateArray[0]->name + ". \nCheck input file");
		}
		// check if the energies and statistical weights are defined
		for (auto& state: stateArray){
			if (state->energy == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " energy for the evaluation of treanorPopulation function. \nCheck input file");
			}
			else if (state->statisticalWeight == Constant::NON_DEF){
				Message::error(std::string("Unable to find ") + state->name + " statistical weight for the evaluation of treanorPopulation function. \nCheck input file");
			}
		}
		// calculate the values of the property
		treanorPopulation(stateArray, argumentArray, argumentStrArray, workCond);
		// assign this function to each state (to be used later on, if necessary)
		for (auto& state: stateArray){
			state->populationFunc = treanorPopulation;
			state->populationParams = argumentArray;
			state->populationStrParams = argumentStrArray;
		}
	}

	else if (propertyString.find("@") != std::string::npos){
		Message::error(std::string("Error! Trying to use the property function '") + propertyFuncStr + "' which is not defined in the code.\n" +
			           "Add it at 'StatePropertyFunctions.h' and do not forget to put it also at the function 'evaluateStatePropertyFunction'.");
	}

	else{
		argumentArray.push_back(Parse::str2value(propertyString));
		constantValue(stateArray, property, argumentArray);
	}
}

#endif
