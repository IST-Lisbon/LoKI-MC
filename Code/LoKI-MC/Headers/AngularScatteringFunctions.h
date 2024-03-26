#ifndef __AngularScatteringFunctions__
#define __AngularScatteringFunctions__
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/MathFunctions.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

// ******************************************************************************  
// The 'AngularScatteringFunctions' generate a cos(random scattering angle), where the scattering angle is between 0 and pi,
// according to the chosen model.                    
// 																				   
// NEW ANGULAR MODELS SHOULD BE ADDED IN THIS FILE!								   
// 																				  
// ******************************************************************************  

namespace AngularScatteringFunctions{
	
	template <class ElectronKineticsType>
	using functionPointer = double (*) (double energy, double energyAfter, int chosenProcessID, ElectronKineticsType* electronKinetics);

	template <class ElectronKineticsType>
	double isotropic(double energy, double energyAfter, int chosenProcessID, ElectronKineticsType* electronKinetics){
		return 1.0-2.0*MathFunctions::unitUniformRand(true,true); // returns cos(angle)
	}

	template <class ElectronKineticsType>
	double forward(double energy, double energyAfter, int chosenProcessID, ElectronKineticsType* electronKinetics){
		return 1; // returns cos(angle)
	}

	template <class ElectronKineticsType>
	double bornDipole(double energy, double energyAfter, int chosenProcessID, ElectronKineticsType* electronKinetics){
		// Based on the Born-Dipole scattering theory presented at
		// L Vialetto 2021: https://doi.org/10.1088/1361-6595/ac0a4d
		// Check equation (25)
		double energyLoss = electronKinetics->energyLosses[chosenProcessID];
		double energyRatio = energyLoss/std::pow(std::sqrt(energyAfter) + std::sqrt(energy), 2);
		double energyRatioSquare = energyRatio*energyRatio;
		return 1.0 + 2.0*energyRatioSquare/(1.0-energyRatioSquare)*(1.0-std::pow(energyRatioSquare, -MathFunctions::unitUniformRand(true,true))); // returns cos(angle)
	}

	template <class ElectronKineticsType>
	double surendra(double energy, double energyAfter, int chosenProcessID, ElectronKineticsType* electronKinetics){
		// Based on the Surendra's scattering theory presented at 
		// Vahedi 1995: https://doi.org/10.1016/0010-4655(94)00171-W 
		// Check equation (9)
		return (2.0 + energy - 2.0*std::pow(1.0+energy,MathFunctions::unitUniformRand(true,true)))/energy; // returns cos(angle)
	}

	template <class ElectronKineticsType>
	double coulombScreen(double energy, double energyAfter, int chosenProcessID, ElectronKineticsType* electronKinetics){
		// Based on Mott's scattering theory presented at
		// Hagelaar 2000: https://iopscience.iop.org/article/10.1088/0963-0252/9/4/318
		// Check equations (21-23)
		std::vector<double> params = electronKinetics->angularScatteringParams[chosenProcessID];
		double epsilon;
		// option 0: consider the incident energy 
		if (params[0] == 0){
			epsilon = energy;
		}
		// option 1: consider the energy after the collision
		else{
			epsilon = energyAfter;
		}
		double screenEnergy = params[1];
		double screenParam = screenEnergy/epsilon;
		double R = MathFunctions::unitUniformRand(true,true);
		return (screenParam + 1.0 - (2.0*screenParam+1.0)*R) / (screenParam + 1.0 - R); // returns cos(angle)
	}

	template <class ElectronKineticsType>
	AngularScatteringFunctions::functionPointer<ElectronKineticsType> functionMap(std::string functionName, std::vector<double> functionParams, ElectronKineticsType* electronKinetics){
		// 'angularScatteringFunctionMap' maps the angular scattering functions
		// NEW ANGULAR SCATTERING FUNCTIONS SHOULD BE ADDED HERE	
		AngularScatteringFunctions::functionPointer<ElectronKineticsType> function;
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
			// not important since this scattering is implemented inside the method 'ionizationCollision' of BoltzmannMC
			function = NULL;
			nParamsRequired = 0;
		}
		else{
			Message::error(std::string("The angular scattering function '") + functionName + "' is not defined in Headers/AngularScatteringFunctions.h.\nDo not forget to add it also at the function 'AngularScatteringFunctionMap'.");
		}
		if (nParamsRequired != functionParams.size()){
			Message::error(std::string("The angular scattering function '") + functionName + "' requires exactly " + std::to_string(nParamsRequired) + " parameters, and not " + std::to_string(functionParams.size()));
		}	
		return function;
	}
};

#endif