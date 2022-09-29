#ifndef __Gas__
#define __Gas__

#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/Constant.h"
#include <iostream>
#include <string>
#include <vector>

class Collision;
class Reaction;
class WorkingConditions;

template <class GasType>
using GasPropertyFunctionPointer = double (*) (GasType* gas, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond);

template <class GasType, class StateType>
class Gas{
	//  Class that stores the information of a certain gas present in the mixture of gases included in the discharge. The electron collisions
	//  and the populations of the different states are to be used in a Boltzmann solver to obtain the EEDF.

public:

	// ----- class attributes -----
	
	int ID=-1; 										// ID that identifies the gas in the gas array
	std::string name; 								// name of the gas
	
	double mass = Constant::NON_DEF; 			    // mass of the gas
	GasPropertyFunctionPointer<GasType> massFunc = NULL;
	std::vector<double> massParams;
	std::vector<std::string> massStrParams;

	double harmonicFrequency = Constant::NON_DEF; 	// harmonic oscillator frequency (molecules)
	GasPropertyFunctionPointer<GasType> harmonicFrequencyFunc = NULL;
	std::vector<double> harmonicFrequencyParams;
	std::vector<std::string> harmonicFrequencyStrParams;

	double anharmonicFrequency = Constant::NON_DEF;	// anharmonic oscillator frequency (molecules)
	GasPropertyFunctionPointer<GasType> anharmonicFrequencyFunc = NULL;
	std::vector<double> anharmonicFrequencyParams;
	std::vector<std::string> anharmonicFrequencyStrParams;

	double rotationalConstant = Constant::NON_DEF;	// rotational constant (molecules)
	GasPropertyFunctionPointer<GasType> rotationalConstantFunc = NULL;
	std::vector<double> rotationalConstantParams;
	std::vector<std::string> rotationalConstantStrParams;

	double lennardJonesDistance = Constant::NON_DEF;// sigma parameter of Lennard-Jones potential
	GasPropertyFunctionPointer<GasType> lennardJonesDistanceFunc = NULL;
	std::vector<double> lennardJonesDistanceParams;
	std::vector<std::string> lennardJonesDistanceStrParams;

	double lennardJonesDepth = Constant::NON_DEF;	// epsilon parameter of Lennard-Jones potential
	GasPropertyFunctionPointer<GasType> lennardJonesDepthFunc = NULL;
	std::vector<double> lennardJonesDepthParams;
	std::vector<std::string> lennardJonesDepthStrParams;

	double electricDipolarMoment = Constant::NON_DEF; // electric dipolar moment (molecules)
	GasPropertyFunctionPointer<GasType> electricDipolarMomentFunc = NULL;
	std::vector<double> electricDipolarMomentParams;
	std::vector<std::string> electricDipolarMomentStrParams;

	double electricQuadrupoleMoment = Constant::NON_DEF; // electric quadrupole moment (molecules)
	GasPropertyFunctionPointer<GasType> electricQuadrupoleMomentFunc = NULL;
	std::vector<double> electricQuadrupoleMomentParams;
	std::vector<std::string> electricQuadrupoleMomentStrParams;

	double polarizability = Constant::NON_DEF;      // polarizability (molecules)
	GasPropertyFunctionPointer<GasType> polarizabilityFunc = NULL;
	std::vector<double> polarizabilityParams;
	std::vector<std::string> polarizabilityStrParams;

	double fraction = 0;                  			// fraction in the gas mixture
	GasPropertyFunctionPointer<GasType> fractionFunc = NULL;
	std::vector<double> fractionParams;
	std::vector<std::string> fractionStrParams;

	double heatCapacity = Constant::NON_DEF;        // constant pressure heat capacity 
	GasPropertyFunctionPointer<GasType> heatCapacityFunc = NULL;
	std::vector<double> heatCapacityParams;
	std::vector<std::string> heatCapacityStrParams;

	double thermalConductivity = Constant::NON_DEF; // fraction in the gas mixture
	GasPropertyFunctionPointer<GasType> thermalConductivityFunc = NULL;
	std::vector<double> thermalConductivityParams;
	std::vector<std::string> thermalConductivityStrParams;

	double OPBParameter = Constant::NON_DEF;
	GasPropertyFunctionPointer<GasType> OPBParameterFunc = NULL;
	std::vector<double> OPBParameterParams;
	std::vector<std::string> OPBParameterStrParams;

	std::vector<StateType*> stateArray;            // State array with all the states of the gas

	// ----- class methods -----

	void disp(){
		std::cout<<"ID: "<<ID<<std::endl;
		std::cout<<"name: "<<name<<std::endl;
		if (mass != Constant::NON_DEF){
			std::cout<<"mass: "<<mass<<std::endl;
		}
		if (harmonicFrequency != Constant::NON_DEF){
			std::cout<<"harmonicFrequency: "<<harmonicFrequency<<std::endl;
		}
		if (anharmonicFrequency != Constant::NON_DEF){
			std::cout<<"anharmonicFrequency: "<<anharmonicFrequency<<std::endl;
		}
		if (rotationalConstant != Constant::NON_DEF){
			std::cout<<"rotationalConstant: "<<rotationalConstant<<std::endl;
		}
		if (lennardJonesDistance != Constant::NON_DEF){
			std::cout<<"lennardJonesDistance: "<<lennardJonesDistance<<std::endl;
		}
		if (lennardJonesDepth != Constant::NON_DEF){
			std::cout<<"lennardJonesDepth: "<<lennardJonesDepth<<std::endl;
		}
		if (electricDipolarMoment != Constant::NON_DEF){
			std::cout<<"electricDipolarMoment: "<<electricDipolarMoment<<std::endl;
		}
		if (electricQuadrupoleMoment != Constant::NON_DEF){
			std::cout<<"electricQuadrupoleMoment: "<<electricQuadrupoleMoment<<std::endl;
		}
		if (polarizability != Constant::NON_DEF){
			std::cout<<"polarizability: "<<polarizability<<std::endl;
		}
		std::cout<<"fraction: "<<fraction<<std::endl;
		if (heatCapacity != Constant::NON_DEF){
			std::cout<<"heatCapacity: "<<heatCapacity<<std::endl;
		}
		if (thermalConductivity != Constant::NON_DEF){
			std::cout<<"thermalConductivity: "<<thermalConductivity<<std::endl;
		}
		if (OPBParameter != Constant::NON_DEF){
			std::cout<<"OPBParameter: "<<OPBParameter<<std::endl;
		}		
	}
	void dispStates(){
		std::cout<<name<<" states:"<<std::endl;
		for(int i=0; i<stateArray.size();++i){
			std::cout<<"\t "<<i<<" - "<<stateArray[i]->name<<std::endl;
		}		
	}
	void dispFamilyTree(){
		std::cout<<name<<" states: (family tree)"<<std::endl;
		std::cout<<name<<std::endl;
		for (auto& eleState: stateArray){
			if (eleState->type == "ele"){
				std::cout<<"\t|"<<eleState->name<<std::endl;
				for(auto& vibState: eleState->childArray){
					std::cout<<"\t|\t|"<<vibState->name<<std::endl;
					for(auto& rotState: vibState->childArray){
						std::cout<<"\t|\t|\t|"<<rotState->name<<std::endl;
					}
				}
			}
		}

		for (auto& ionState: stateArray){
			if (ionState->type=="ion"){
				std::cout<<"\t|"<<ionState->name<<std::endl;
			}
		}
	}
	void evaluateDensities(){
	    // evaluateDensities is a function that evaluates the densities of all the states of the gas from their
	    // populations and the gas fraction. Note that the densities of the states are normalized to the total gas
	    // density.
	    for (int i=0; i<stateArray.size();++i){
	    	stateArray[i]->evaluateDensity();
	    }
	} 
	static int add(std::string gasName, std::vector<GasType*> &gasArray){
	    // addGas looks for gasName in Gas array gas. If a Gas with the same
	    // name is found the function returns its ID. If no Gas is found with
	    // the same name, a new Gas is added to the array and its ID is returned.

	    int gasID = Gas::find(gasName, gasArray);
	    if (gasID != -1){
	    	return gasID;
	    }
	    else{
	    	gasArray.push_back(  new GasType(gasName));
	    	return gasArray.back()->ID;
	    }		
	}
	static int find(std::string gasName,std::vector<GasType*> gasArray){

		for (auto& gas: gasArray){
			if (gasName == gas->name){
				return gas->ID;
			}
		}
		return -1;
	}
	static void checkFractionNorm(std::vector<GasType*> gasArray){
	    // checkFractionsNorms checks for the fractions of the different gases in gasArray to add to one. In case the fractions of the gases are not
	    // properly normalised, an error message is thrown.

	    double norm=0;
	    for (auto& gas: gasArray){
	    	norm = norm + gas->fraction;
	    }
	    if (std::abs(norm-1)>10*std::numeric_limits<double>::epsilon()){
	    	Message::error(std::string("Gas fractions are not properly normalized (Error = ") + std::to_string(norm-1) + "). Please, check input file.\n");
	    }		
	}
};

#endif