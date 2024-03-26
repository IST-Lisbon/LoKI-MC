#ifndef __AngularDistributionFunctions__
#define __AngularDistributionFunctions__
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Collision.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

// ******************************************************************************  
// The 'AngularDistributionFunctions', I(epsilon,theta), represent the probability 
// of having a scattered electron with a certain angle theta, when the incident    
// electron has an energy epsilon. These functions are normalized such that:       
// 2*Pi*Integral_0^Pi (I(epsilon,theta)*sin(theta) dTheta) = 1                     
// 																				   
// NEW ANGULAR MODELS SHOULD BE ADDED IN THIS FILE!								   
// 																				  
// These distribution functions are used to obtain the momentum-transfer cross sections
// (check documentation).														
// ******************************************************************************  						


namespace AngularDistributionFunctions{
	
	/*
		'funcionPointer' is defined in Collision.h in the following way:
		using functionPointer = double (*) (double energy, double angle, Collision* collision, bool superElasticComponent);
	*/

	static double isotropic(double energy, double angle, Collision* collision, bool superElasticComponent){
		// not important, since it is not used in the code
		// When using isotropic, momTransCrossSection = integralCrossSection
		return 1.0/(4.0*M_PI);
	}

	static double forward(double energy, double angle, Collision* collision, bool superElasticComponent){
		// not important, since it is not used in the code
		// When using forward, momTransCrossSection = 0
		double value;
		if (angle == 0){
			value = 1.0;
		}
		else{
			value = 0;
		}
		return value;
	}

	static double bornDipole(double energy, double angle, Collision* collision, bool superElasticComponent){
		// Based on the Born-Dipole scattering theory presented at
		// L Vialetto 2021: https://doi.org/10.1088/1361-6595/ac0a4d
		// Check equation (23)
		double energyAfter;
		if (superElasticComponent){
			energyAfter = energy + collision->threshold;
		} 
		else{
			energyAfter = energy - collision->threshold;
		}
		if (energyAfter <= 0){
			return 1.0/(4.0*M_PI);
		}	
		double momentumChangeSquared = energyAfter + energy - 2.0*std::sqrt(energy*energyAfter)*std::cos(angle);
		double sqrtEnergy = std::sqrt(energy), sqrtEnergyAfter = std::sqrt(energyAfter);
		return 1.0/(4.0*M_PI)*sqrtEnergy*sqrtEnergyAfter/momentumChangeSquared / 
			std::log((sqrtEnergyAfter+sqrtEnergy)/std::sqrt(collision->threshold));
	}

	static double surendra(double energy, double angle, Collision* collision, bool superElasticComponent){
		// Based on the Surendra's scattering theory presented at 
		// Vahedi 1995: https://doi.org/10.1016/0010-4655(94)00171-W 
		// Check equation (7)
		if (energy == 0){
			return 1.0/(4.0*M_PI);
		}
		return energy / (4.0*M_PI *(1.0+energy*std::pow(std::sin(angle/2.0),2)) * std::log(1.0+energy));
	}

	static double coulombScreen(double energy, double angle, Collision* collision, bool superElasticComponent){
		// Based on Mott's scattering theory presented at
		// Hagelaar 2000: https://iopscience.iop.org/article/10.1088/0963-0252/9/4/318
		// Check equations (21-23)
		double epsilon;
		// option 0: consider the incident energy 
		if (collision->angularScatteringParams[0] == 0){
			epsilon = energy;
		}
		// option 1: consider the energy after the collision 
		else{
			if (superElasticComponent){
				epsilon = energy + collision->threshold;
			} 
			else{
				epsilon = energy - collision->threshold;
			}
		}
		if (epsilon <= 0){
			return 1.0/(4.0*M_PI);
		}	
		double screenEnergy = collision->angularScatteringParams[1];
		double screenParam = screenEnergy/epsilon;
		return (screenParam*(screenParam+1.0))/M_PI/std::pow(2.0*screenParam+1.0-std::cos(angle),2);
	}

	static double momentumConservationIonization(double energy, double angle, Collision* collision, bool superElasticComponent){
		// not true (of course). TO BE IMPROVED
		// it depends also on the type of energy-sharing!
		// This does not affect the calculation of the electron kinetics, only the total MT cross-section used for the calculation of swarm params from the EEDF
		return 1.0/(4.0*M_PI);
	}

	static functionPointer functionMap(std::string functionName, std::vector<double> functionParams){
		// 'FunctionMap' maps the angular distribution functions
		// NEW ANGULAR DISTRIBUTION FUNCTIONS SHOULD BE ADDED HERE
		functionPointer function;
		int nParamsRequired = -1;
		if (functionName == "isotropic"){
			function = isotropic;
			nParamsRequired = 0;
		}
		else if (functionName == "forward"){
			function = forward;
			nParamsRequired = 0;
		}
		else if (functionName == "bornDipole"){
			function = bornDipole;
			nParamsRequired = 0;
		}
		else if (functionName == "surendra"){
			function = surendra;
			nParamsRequired = 0;
		}
		else if (functionName == "coulombScreen"){
			function = coulombScreen;
			nParamsRequired = 2;
		}
		else if (functionName == "momentumConservationIonization"){
			function = momentumConservationIonization;
			nParamsRequired = 0;
		}
		else{
			Message::error(std::string("The angular distribution function '") + functionName + "' is not defined in Headers/AngularDistributionFunctions.h.\nDo not forget to add it also at the function 'angularDistributionFunctionMap'.");
		}
		if (nParamsRequired != functionParams.size()){
			Message::error(std::string("The angular distribution function '") + functionName + "' requires exactly " + std::to_string(nParamsRequired) + " parameters, and not " + std::to_string(functionParams.size()));
		}
		return function;
	}
};

#endif